/********************  HEADERS  *********************/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include "lbm_comm.h"


/*******************  FUNCTION  *********************/
/**
 * Affiche la configuation du lbm_comm pour un rank donné
 * @param mesh_comm Configuration à afficher
**/
void  lbm_comm_print( lbm_comm_t *mesh_comm )
{
	int rank ;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	printf( " RANK %d (TOP %d BOTTOM %d) ( POSITION %d %d ) (WH %d %d ) \n", rank,
										mesh_comm->top_id,
									    mesh_comm->bottom_id,
									    mesh_comm->x,
									    mesh_comm->y,
									    mesh_comm->width,
									    mesh_comm->height );
}


/*******************  FUNCTION  *********************/
/**
 * Initialise un lbm_comm :
 * - Voisins
 * - Taille du maillage local
 * - Position relative
 * @param mesh_comm MeshComm à initialiser
 * @param rank Rank demandant l'initalisation
 * @param comm_size Taille totale du communicateur
 * @param width largeur du maillage
 * @param height hauteur du maillage
**/
void lbm_comm_init( lbm_comm_t * mesh_comm, int rank, int comm_size, int width, int height )
{
	//vars
	int nb_x;
	int nb_y;
	int nlines;
	int extras;

	//compute splitting
	nb_y = comm_size <= height ? comm_size : 0;
	nb_x = 1;

	if (nb_y == 0)
		fatal("Can't get a 2D cut for current problem size and number of processes.");

	nlines = height / comm_size;
	extras = height % comm_size;


	//setup nb
	mesh_comm->nb_x = nb_x;
	mesh_comm->nb_y = nb_y;

	//setup size (+2 for ghost cells on border)
	mesh_comm->width  = width+2;
	mesh_comm->height = rank == 0 ? nlines+extras+2 : nlines+2;

	//setup position
	mesh_comm->x = 0;
	mesh_comm->y = rank * nlines + extras;
	
	// Compute neighbour nodes id
	mesh_comm->top_id     = rank == 0           ? -1 : rank-1;
	mesh_comm->bottom_id  = rank == comm_size-1 ? -1 : rank+1;

	//if debug print comm
	#ifndef NDEBUG
	lbm_comm_print(mesh_comm);
	#endif
}


/*******************  FUNCTION  *********************/
/**
 * Libere un lbm_comm
 * @param mesh_comm MeshComm à liberer
**/
void lbm_comm_release( lbm_comm_t * mesh_comm )
{
	mesh_comm->x = 0;
	mesh_comm->y = 0;
	mesh_comm->width = 0;
	mesh_comm->height = 0;
	mesh_comm->top_id = -1;
	mesh_comm->bottom_id = -1;
}


void lbm_comm_ghost_send(lbm_comm_t *mesh, Mesh *mesh_to_process)
{
	MPI_Request req[2];
	int size = mesh_to_process->width*DIRECTIONS;

	if(mesh->bottom_id != -1) {
		MPI_Isend(Mesh_get_cell(mesh_to_process, 1, mesh->height-2), size, MPI_DOUBLE, mesh->bottom_id, 0, MPI_COMM_WORLD, req);
		//MPI_Wait(req, MPI_STATUS_IGNORE);
	}

	if(mesh->top_id != -1) {
		MPI_Isend(Mesh_get_cell(mesh_to_process, 1, 1), size, MPI_DOUBLE, mesh->top_id, 0, MPI_COMM_WORLD, req+1);
		//MPI_Wait(req+1, MPI_STATUS_IGNORE);
	}
}

void lbm_comm_ghost_recv(lbm_comm_t *mesh, Mesh *mesh_to_process)
{
	MPI_Status status;
	int size = mesh_to_process->width*DIRECTIONS;

	if(mesh->bottom_id != -1)
		MPI_Recv(Mesh_get_cell(mesh_to_process, 1, mesh->height-1), size, MPI_DOUBLE, mesh->bottom_id, 0, MPI_COMM_WORLD, &status);

	if(mesh->top_id != -1)
		MPI_Recv(Mesh_get_cell(mesh_to_process, 1, 0), size, MPI_DOUBLE, mesh->top_id, 0, MPI_COMM_WORLD, &status);
}

void lbm_comm_ghost_exchange(lbm_comm_t *mesh, Mesh *mesh_to_process)
{
	MPI_Status status;
	MPI_Request req[2];
	int size = mesh_to_process->width*DIRECTIONS;

	if(mesh->bottom_id != -1) {
		MPI_Isend(Mesh_get_cell(mesh_to_process, 1, mesh->height-2), size, MPI_DOUBLE, mesh->bottom_id, 0, MPI_COMM_WORLD, req);
		MPI_Wait(req, MPI_STATUS_IGNORE);
	}

	if(mesh->top_id != -1) {
		MPI_Isend(Mesh_get_cell(mesh_to_process, 1, 1), size, MPI_DOUBLE, mesh->top_id, 0, MPI_COMM_WORLD, req+1);
		MPI_Wait(req+1, MPI_STATUS_IGNORE);
	}


	if(mesh->bottom_id != -1)
		MPI_Recv(Mesh_get_cell(mesh_to_process, 1, mesh->height-1), size, MPI_DOUBLE, mesh->bottom_id, 0, MPI_COMM_WORLD, &status);

	if(mesh->top_id != -1)
		MPI_Recv(Mesh_get_cell(mesh_to_process, 1, 0), size, MPI_DOUBLE, mesh->top_id, 0, MPI_COMM_WORLD, &status);
}

/*******************  FUNCTION  *********************/
/**
 * Rendu du mesh en effectuant une réduction a 0
 * @param mesh_comm MeshComm à utiliser
 * @param temp Mesh a utiliser pour stocker les segments
**/
void save_frame_all_domain( FILE * fp, Mesh *source_mesh, Mesh *temp )
{
	//vars
	int i = 0;
	int comm_size, rank ;
	MPI_Status status;

	//get rank and comm size
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	/* If whe have more than one process */
	if(1 < comm_size)
	{
		if(rank == 0)
		{
			/* Rank 0 renders its local Mesh */
			save_frame(fp,source_mesh);
			/* Rank 0 receives & render other processes meshes */
			for(i = 1 ; i < comm_size ; i++)
			{
				MPI_Recv(temp->cells, source_mesh->width*source_mesh->height*DIRECTIONS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				save_frame(fp,temp);
			}
		} else {
			/* All other ranks send their local mesh */
			MPI_Send( source_mesh->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		}
	} else {
		/* Only 0 renders its local mesh */
		save_frame(fp,source_mesh);
	}

}


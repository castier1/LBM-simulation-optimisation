/********************  HEADERS  *********************/
#include <mpi.h>
#include <omp.h>
#include <assert.h>
#include "lbm_phys.h"
#include "lbm_init.h"

/*******************  FUNCTION  *********************/
/**
 * Applique les conditions initiale avec un fluide au repos. Utiliser uniquement
 * pour des tests lors de l'implémentation.
**/
//Apply initial conditions :
// - v = 0
// - rho = 1
//Equation is : f_equilibre(x,t) = w_i * rho + rho * s_i(v(x,t))
//Here v = 0 so s_i(*) = 0, rho = 1, so keep only w_i
//w_i is the direction weight.
void init_cond_velocity_0_density_1(Mesh * mesh)
{
	//vars
	int i,j,k;

	//errors
	assert(mesh != NULL);

#pragma omp for private(j,k)
	//loop on all cells
	for ( i = 0 ; i <  mesh->height ; i++)
		for ( j = 0 ; j <  mesh->width ; j++)
			for ( k = 0 ; k < DIRECTIONS ; k++)
				Mesh_get_cell(mesh, j, i)[k] = equil_weight[k];
}

/*******************  FUNCTION  *********************/
/**
 * Initialisation de l'obstacle, on bascule les types des mailles associés à CELL_BOUNCE_BACK.
 * Ici l'obstacle est un cercle de centre (OBSTACLE_X,OBSTACLE_Y) et de rayon OBSTACLE_R.
**/
void setup_init_state_circle_obstacle(Mesh * mesh, lbm_mesh_type_t * mesh_type, const lbm_comm_t * mesh_comm)
{
	//vars
	int i,j;

	//loop on nodes
#pragma omp for private(j)
	for (i = mesh_comm->y; i < mesh->height + mesh_comm->y; i++)
		for (j = mesh_comm->x; j < mesh->width + mesh_comm->x; j++)
			if (((i-OBSTACLE_Y) * (i-OBSTACLE_Y)) + ((j-OBSTACLE_X) * (j-OBSTACLE_X)) <= OBSTACLE_R*OBSTACLE_R)
				*(lbm_cell_type_t_get_cell(mesh_type, j-mesh_comm->x, i-mesh_comm->y)) = CELL_BOUNCE_BACK;
}

/*******************  FUNCTION  *********************/
/**
 * Initialise le fluide complet avec un distribution de poiseuille correspondant un état d'écoulement
 * linéaire à l'équilibre.
 * @param mesh Le maillage à initialiser.
 * @param mesh_type La grille d'information notifiant le type des mailles.
**/
void setup_init_state_global_poiseuille_profile(Mesh * mesh, lbm_mesh_type_t * mesh_type,const lbm_comm_t * mesh_comm)
{
	//vars
	int i,j,k;
	Vector v = {0.0,0.0};
	const double density = 1.0;

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//apply poiseuil for all nodes except on top/bottom border
#pragma omp for private(j,k)
	for ( i = 0 ; i < mesh->height ; i++)
	{
		for ( j = 0 ; j < mesh->width ; j++)
		{
			for ( k = 0 ; k < DIRECTIONS ; k++)
			{
				//compute equilibr.
				v[0] = helper_compute_poiseuille(i + mesh_comm->y,MESH_HEIGHT);
				Mesh_get_cell(mesh, j, i)[k] = compute_equilibrium_profile(v,density,k);
				//mark as standard fluid
				*( lbm_cell_type_t_get_cell( mesh_type , j, i) ) = CELL_FUILD;
				//this is a try to init the fluide with null speed except on left interface.
				if (j > 1)
					Mesh_get_cell(mesh, j, i)[k] = equil_weight[k];
			}
		}
	}
}

/*******************  FUNCTION  *********************/
/**
 * Initialisation des conditions aux bords.
 * @param mesh Le maillage à initialiser.
 * @param mesh_type La grille d'information notifiant le type des mailles.
**/
void setup_init_state_border(Mesh * mesh, lbm_mesh_type_t * mesh_type, const lbm_comm_t * mesh_comm)
{
	//vars
	int i,k;
	Vector v = {0.0,0.0};
	const double density = 1.0;

	//top
	if (mesh_comm->top_id == -1)
#pragma omp for private(k)
		for ( i = 0 ; i < mesh->width ; i++)
			for ( k = 0 ; k < DIRECTIONS ; k++)
			{
				//compute equilibr.
				Mesh_get_cell(mesh, i, 0)[k] = compute_equilibrium_profile(v,density,k);
				//mark as bounce back
				*( lbm_cell_type_t_get_cell(mesh_type, i, 0) ) = CELL_BOUNCE_BACK;
			}

	//bottom
	if (mesh_comm->bottom_id == -1)
#pragma omp for private(k)
		for ( i = 0 ; i < mesh->width ; i++)
			for ( k = 0 ; k < DIRECTIONS ; k++)
			{
				//compute equilibr.
				Mesh_get_cell(mesh, i, mesh->height - 1)[k] = compute_equilibrium_profile(v,density,k);
				//mark as bounce back
				*( lbm_cell_type_t_get_cell(mesh_type , i, mesh->height - 1)) = CELL_BOUNCE_BACK;
			}
}

/*******************  FUNCTION  *********************/
/**
 * Mise en place des conditions initiales.
 * @param mesh Le maillage à initialiser.
 * @param mesh_type La grille d'information notifiant le type des mailles.
 * @param mesh_comm La structure de communication pour connaitre notre position absolue dans le maillage globale.
**/
void setup_init_state(Mesh * mesh, lbm_mesh_type_t * mesh_type, const lbm_comm_t * mesh_comm)
{
#pragma omp parallel
	{
		setup_init_state_global_poiseuille_profile(mesh,mesh_type,mesh_comm);
		setup_init_state_border(mesh,mesh_type,mesh_comm);
		setup_init_state_circle_obstacle(mesh,mesh_type, mesh_comm);
	}
}

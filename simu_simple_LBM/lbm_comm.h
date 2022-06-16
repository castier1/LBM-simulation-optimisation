#ifndef MESH_COMM_H
#define MESH_COMM_H

/********************  HEADERS  *********************/
#include <mpi.h>
#include <stdlib.h>
#include "lbm_struct.h"

/*******************  DEFINITIONS  ******************/
/** Definition de l'ID du processus maître. **/
#define RANK_MASTER 0

/********************  STRUCT  **********************/
/**
 * Structure utilisée pour stoquer les informations relatives aux communications.
**/
typedef struct lbm_comm_t_s
{
	/** Position de la maille locale dans le maillage global (origine). **/
	int x;
	int y;
	/** Taille de la maille locale. **/
	int width;
	int height;
	int nb_x;
	int nb_y;
	int top_id;
	int bottom_id;
} lbm_comm_t;

/*******************  FUNCTION  *********************/
static inline int lbm_comm_width( lbm_comm_t *mc )
{
	return mc->width;
}

/*******************  FUNCTION  *********************/
static inline int lbm_comm_height( lbm_comm_t *mc )
{
	return mc->height;
}

/*******************  FUNCTION  *********************/
void lbm_comm_init( lbm_comm_t * mesh, int rank, int comm_size, int width, int height );
void lbm_comm_release( lbm_comm_t * mesh );
void lbm_comm_print( lbm_comm_t *mesh );

/*******************  FUNCTION  *********************/
void lbm_comm_ghost_exchange(lbm_comm_t * mesh, Mesh *mesh_to_process );
void lbm_comm_ghost_send(lbm_comm_t *mesh, Mesh *mesh_to_process);
void lbm_comm_ghost_recv(lbm_comm_t *mesh, Mesh *mesh_to_process);

/*******************  FUNCTION  *********************/
void save_frame_all_domain( FILE * fp, Mesh *source_mesh, Mesh *temp );

#endif

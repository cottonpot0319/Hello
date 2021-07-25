/******************************************************************************
*         Allocate function                                                   *
*******************************************************************************/
extern void* allocate_vector  (int, int);
extern void  deallocate_vector(void*);
extern void** allocate_matrix (int, int, int);
extern void deallocate_matrix (void**);
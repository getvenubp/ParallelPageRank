/****************************************************************************

 Object version :
    Sang Lim(slim@npac.syr.edu)
    Northeast Parallel Architectures Center at Syracuse University
    11/28/98
****************************************************************************/

import mpi.*;

class alltoallO {
  static public void main(String[] args) throws MPIException {

    final int MAXLEN = 10;

    int i,j,k,l;
    int myself,tasks;
 
    MPI.Init(args);
    myself = MPI.COMM_WORLD.Rank();
    tasks = MPI.COMM_WORLD.Size();

    int out[][] = new int[MAXLEN*tasks][MAXLEN]; 
    int in[][]  = new int[MAXLEN*tasks][];

    for(k=0;k<MAXLEN;k++)
      for(i=0;i<MAXLEN*tasks;i++)  
        out[i][k] = k;   
    
    MPI.COMM_WORLD.Alltoall(out,0,MAXLEN,MPI.OBJECT,in,0,MAXLEN,MPI.OBJECT);

    for(k=0;k<MAXLEN*tasks;k++) 
      for(l=0; l<MAXLEN; l++){
 	if(in[k][l] != l) {
	  System.out.println("bad answer in["+k+"]["+l+"] = "+
		  	       in[k][l]+" should be "+l);
          break; 
	}
      }

    MPI.COMM_WORLD.Barrier();
    if(myself==0)  System.out.println("AllToAllO TEST COMPLETE\n");
    MPI.Finalize();
  }
}



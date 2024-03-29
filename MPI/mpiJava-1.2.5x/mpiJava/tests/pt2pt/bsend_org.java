/****************************************************************************

 MESSAGE PASSING INTERFACE TEST CASE SUITE

 Copyright IBM Corp. 1995

 IBM Corp. hereby grants a non-exclusive license to use, copy, modify, and
 distribute this software for any purpose and without fee provided that the
 above copyright notice and the following paragraphs appear in all copies.

 IBM Corp. makes no representation that the test cases comprising this
 suite are correct or are an accurate representation of any standard.

 In no event shall IBM be liable to any party for direct, indirect, special
 incidental, or consequential damage arising out of the use of this software
 even if IBM Corp. has been advised of the possibility of such damage.

 IBM CORP. SPECIFICALLY DISCLAIMS ANY WARRANTIES INCLUDING, BUT NOT LIMITED
 TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS AND IBM
 CORP. HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 ENHANCEMENTS, OR MODIFICATIONS.

****************************************************************************

 These test cases reflect an interpretation of the MPI Standard.  They are
 are, in most cases, unit tests of specific MPI behaviors.  If a user of any
 test case from this set believes that the MPI Standard requires behavior
 different than that implied by the test case we would appreciate feedback.

 Comments may be sent to:
    Richard Treumann
    treumann@kgn.ibm.com

****************************************************************************

 MPI-Java version :
    Sung-Hoon Ko(shko@npac.syr.edu)
    Northeast Parallel Architectures Center at Syracuse University
    03/22/98

****************************************************************************
*/

import mpi.*;

class bsend {
  static public void main(String[] args) throws MPIException {
    /* Note that the buffer sizes must include the BSEND_OVERHEAD;
       these values are probably sizeof(int) too large */
    int len,tsks,me,i,size,rc;	
    final int A1   = 1000;
    final int A100 = 100000;
    Status status;

    MPI.Init(args);
    me=MPI.COMM_WORLD.Rank();


    int data1[] = new int[A1];
    int data100[] = new int[A100];
    //int intsize = 4;
    int intsize = MPI.INT.Size();

    byte buf1[] = new byte[A1*intsize+MPI.BSEND_OVERHEAD];
    byte buf100[] = new byte[A100*intsize+MPI.BSEND_OVERHEAD];

    if ( me == 0 ) {      
      MPI.Buffer_attach(buf1);  
      MPI.COMM_WORLD.Bsend(data1,0,A1,MPI.INT,1,1);
      MPI.Buffer_detach(buf1);
      
      MPI.Buffer_attach(buf100);
      
      MPI.COMM_WORLD.Barrier();
      
      /* test to see if large array is REALLY being buffered */
      for(i=0;i<A100;i++)  data100[i] = 1;
      
      MPI.COMM_WORLD.Bsend(data100,0,A100,MPI.INT,1,1);
      
      MPI.COMM_WORLD.Recv(data100,0,A100,MPI.INT,1,2);
      
      for(i=0;i<A100;i++)
	if(data100[i] != 2)  
	  System.out.println
	    ("ERROR, incorrect data["+i+"]="+data100[i]+", task 0");
    } else if ( me == 1 ) {
      MPI.COMM_WORLD.Recv(data1,0,A1,MPI.INT,0,1);
      MPI.COMM_WORLD.Barrier();
      
      MPI.Buffer_attach(buf100);
      
      /* test to see if large array is REALLY being buffered */
      for(i=0;i<A100;i++)  data100[i] = 2;
      
      MPI.COMM_WORLD.Bsend(data100,0,A100,MPI.INT,0,2);
      
      MPI.COMM_WORLD.Recv(data100,0,A100,MPI.INT,0,1);
      
      for(i=0;i<A100;i++)
	if(data100[i] != 1)  
	  System.out.println
	    ("ERROR , incorrect data["+i+"]="+data100[i]+", task 1");
      
    }

    MPI.COMM_WORLD.Barrier();
    if(me == 0)  System.out.println("Bsend TEST COMPLETE\n");
    MPI.Finalize();
  }

}

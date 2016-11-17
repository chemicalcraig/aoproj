#ifndef UTIL_H
#define UTIL_H

#define ang2au 1.889725989;
/** Kronecker Delta **/
#define Kronecker(a,b) ((a == b ? 1 : 0))

#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include "atom.h"
#include "molecule.h"
#include <iomanip>
#include <omp.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
using namespace std;

/* Search strings */
string geom_str=" Output coordinates";
string basislabel_str="    Basis function labels";
string overlap_str = " global array: AO ovl";
string nbasis_str="          This is a Direct SCF calculation.";
string tddft_str="  Convergence criterion met";
string tddft_state_stop_str="-----------------------------";
string tddft_task_str1="tddft";
string tddft_task_str2="TDDFT";
string tddft_task_str3="Tddft";
string mo_analysis_str="                       DFT Final Molecular Orbital Analysis";
string lindep_str = " !! The overlap matrix has";
string aobas_str = "          AO basis - number of functions:";

/* Number of columns for printing of matrices */
int ncol = 6;

/********************************************************
 * Functions
 * *******************************************************/
double* d_app_array(double *a,double x) {
  int s1=sizeof(a)/sizeof(double);
  double *temp = new double[s1+1];
  for (int i=0; i<s1; i++) {
    temp[i] = a[i];
  }
  temp[s1] = x;
  return temp;
}

/** Calculate WF overlap **/
bool calcOverlap2(Molecule *mol, int root,Molecule *mol2, int root2) {
  double dum=0.;
//  for (int i=0; i<mol->roots[root].nci; i++) {
//    for (int j=0; j<mol2->roots[root].nci; j++) {
//      dum +=  mol->roots[root].ci[i]
//              *mol2->roots[root2].ci[j]
//              *Kronecker(i,j);
//    }
//  }
//  cout<<"Maximal overlap = "<<dum<<endl;

  /* Using truncated CI expansions, threshold set by NWChem */
  /* for reference, see JCTC, 12, 1207, 2016 */
  /* for open shell, we're only doing alpha electronic part,
   * as the alpha and beta orbitals are the same, and the 
   * alpha and beta CI coefficients are the same for s^2=0.
   * for S^2 = 2...
   * The overlap is multiplied by 2 to get the correct answer*/
  double sum = 0.; 
  double s2,s3,s4;
  s2=0.;
  s3=0.;
  s4=0.;
  int active=20;
  double *alphamatrix = new double[mol->nbasis*mol2->nbasis];
  double thresh=1.e-5;
  
  /* Get number of non-zero CI coefficients for each molecule */
  int nci1=0;
  int nci2=0;
  int nci=0;
  for (int i=0; i<mol->nocc; i++) {
    for (int j=0; j<mol->nbasis; j++) {
      if (fabs(mol->roots[root].ci[nci1]) >= thresh) {
        nci1 ++;
      }
      if (fabs(mol2->roots[root2].ci[nci2]) >= thresh) {
        nci2 ++;
      }
    }
  }
  nci = 2*(nci1+nci2);
  /** Save the non-zero ci coefficient indices **/
  int ci1[nci][2],ci2[nci][2];
  nci1=0;
  nci2=0;

  if (mol->os) {
  cout<<"open shell"<<endl;
  for (int i=0; i<mol->nocc; i++) {
    for (int j=mol->nocc; j<mol->nbasis; j++) {
      if (fabs(mol->cia[i+j*mol->nmo+root*mol->nmo*mol->nmo]) >= thresh) {
        ci1[nci1][0] = i;
        ci1[nci1][1] = j;
        //cout<<"NCI1 "<<nci1<<" "<<i<<" "<<j<<endl;
        nci1 ++;
      }
      if (fabs(mol2->cia[i+j*mol2->nmo+root2*mol2->nmo*mol2->nmo]) >= thresh) {
        ci2[nci2][0] = i;
        ci2[nci2][1] = j;
       // cout<<"NCI2 "<<nci2<<" "<<i<<" "<<j<<endl;
        nci2 ++;
      }
    }
  } 
  for (int i=0; i<mol->nocc; i++) {
    for (int j=mol->nocc; j<mol->nbasis; j++) {
      if (fabs(mol->cib[i+j*mol->nmo+root*mol->nmo*mol->nmo]) >= thresh) {
        ci1[nci1][0] = i;
        ci1[nci1][1] = j;
     //   cout<<"NCI1 "<<nci1<<" "<<i<<" "<<j<<endl;
        nci1 ++;
      }
      if (fabs(mol2->cib[i+j*mol2->nmo+root2*mol2->nmo*mol2->nmo]) >= thresh) {
        ci2[nci2][0] = i;
        ci2[nci2][1] = j;
       // cout<<"NCI2 "<<nci2<<" "<<i<<" "<<j<<endl;
        nci2 ++;
      }
    }
  } 
  }//end os
  else {
  for (int i=0; i<mol->nocc; i++) {
    for (int j=mol->nocc; j<mol->nbasis; j++) {
      if (fabs(mol->ci[i+j*mol->nmo+root*mol->nmo*mol->nmo]) >= thresh) {
        ci1[nci1][0] = i;
        ci1[nci1][1] = j;
     //   cout<<"NCI1 "<<nci1<<" "<<i<<" "<<j<<endl;
        nci1 ++;
      }
      if (fabs(mol2->ci[i+j*mol2->nmo+root2*mol2->nmo*mol2->nmo]) >= thresh) {
        ci2[nci2][0] = i;
        ci2[nci2][1] = j;
       // cout<<"NCI2 "<<nci2<<" "<<i<<" "<<j<<endl;
        nci2 ++;
      }
    }
  } 
  } //end closed shell

  int s;
  double det,det2;
  gsl_matrix *mat = gsl_matrix_calloc(active,active);
  gsl_matrix *mat2 = gsl_matrix_calloc(active,active);
  gsl_matrix *mol1mo = gsl_matrix_calloc(mol->nmo,mol->nmo);
  gsl_matrix *mol2mo = gsl_matrix_calloc(mol->nmo,mol->nmo);
  gsl_matrix *ludecomp = gsl_matrix_calloc(active,active);
  gsl_permutation *perm = gsl_permutation_calloc(mat->size1);
  
  bool first=true;
  int from1,from2,to1,to2;
  int mo1,mo2;

  for (int I=0; I<mol->roots[root].ncia; I++) { //CI 1
    for (int J=0; J<mol2->roots[root2].ncia; J++) { //CI 2
      from1 = ci1[I][0];
      to1 = ci1[I][1];
      from2 = ci2[J][0];
      to2 = ci2[J][1];
      sum=0.;
      
      /** copy mo matrices into gsl format **/
      /** Assumes mol1 and mol2 have equal nmo and nbasis **/
      for (int i=0; i<mol->nmo; i++)
        for (int j=0; j<mol->nbasis; j++) {
          gsl_matrix_set(mol1mo,i,j,mol->mos[i+j*mol->nmo]);
          gsl_matrix_set(mol2mo,i,j,mol2->mos[i+j*mol2->nmo]);
          //cout<<i<<" "<<j<<" "<<gsl_matrix_get(mol1mo,i,j)<<endl;
        }
      
      /** Swap columns in MO matrices **/
      cout<<from1<<" "<<to1<<" "<<from2<<" "<<to2<<endl;
      //exit(0);
      gsl_matrix_swap_rows(mol1mo,from1,to1);
      gsl_matrix_swap_rows(mol2mo,from2,to2);
    
    double sum2 = 0.;
    /** Build alpha matrix and get determinant **/
      for (int i=mol->nocc-active; i<mol->nocc/*2*active*/; i++) { //MO 1
        for (int j=mol2->nocc-active; j<mol2->nocc/*2*active*/; j++) { //MO 2
          mo1 = i;
          mo2 = j;
          sum=0.;
          sum2 = 0.;

         #pragma omp parallel for reduction (+:sum,sum2)
          for (int a=0; a<mol->nbasis; a++) { //AO 1
            for (int b=0; b<mol2->nbasis; b++) { //AO 2    
              /** Generate alpha matrix **/
              sum += gsl_matrix_get(mol1mo,mo1,a)
                * gsl_matrix_get(mol2mo,mo2,b)
                * mol->overlapm[a+b*mol->nbasis];
              
              sum2 += mol->mos[mo1+a*mol->nmo]
                  * mol2->mos[mo2+b*mol2->nmo]
                  * mol->overlapm[a+b*mol->nbasis];
            
            } //AO2
          } //AO1

          int row = i-(mol->nocc-active);
          int col = j-(mol->nocc-active);
          gsl_matrix_set(mat,row,col,sum);
          gsl_matrix_set(mat2,row,col,sum2);
          
          double sum3 = 1.;
          for (int x=0; x<active; x++) {
            sum3 *= gsl_matrix_get(mat,x,x);
            //cout<<"diag "<<x<<" "<<gsl_matrix_get(mat,x,x)<<" "<<sum3<<endl;
          }
          //cout<<"diag "<<sum3<<endl;
          //cout<<"gsl "<<i<<" "<<j<<" "<<gsl_matrix_get(mat,row,col)<<" "<<gsl_matrix_get(mat2,row,col)<<endl;
        }//MO2
      
      //exit(0);
      }//MO1
      /** Unperturbed Alpha matrix complete **/
      gsl_matrix_memcpy(ludecomp, mat);
      /** take determinant of alpha matrix **/
      gsl_linalg_LU_decomp(ludecomp, perm, &s);
      det=gsl_linalg_LU_det(ludecomp,s);
      if (first) {
        gsl_matrix_memcpy(ludecomp, mat2);
        /** take determinant of alpha matrix **/
        gsl_linalg_LU_decomp(ludecomp, perm, &s);
        det2=gsl_linalg_LU_det(ludecomp,s);

        first = false;
      } //end if first
      
      double d1 = mol->roots[root].ci[I];
      double d2 = mol2->roots[root2].ci[J];
      double d3 =  d1*d2;

      s2 += d3;
      s3 += d1*d2*det*det2;
      cout<<"("<<I<<","<<J<<"): |alpha|="<<det<<",    |beta|="<<det2<<",  c_I,J = "<<d1<<", "<<d2
      <<",    s(I,J) = "<<s3*2<<endl;
    }//CI2
  }//CI1
  mol->projection = s3*2;

  sum=0.;
  for (int i=0; i<mol->roots[root].nci; i++) {
    sum += mol->roots[root].ci[i]*mol->roots[root].ci[i];
  }
  cout<<"S("<<root+1<<","<<root2+1<<") is "<<100*mol->projection/sum<<"% of maximal overlap"<<endl;
  return true;
}


/**** Skip several lines in input file ****/
void getnlines(ifstream &in, char *temp, int n, int length) {
  for (int i=0; i<n; i++) {
    in.getline(temp,length);
  }
}

bool calcdip(Molecule *mol) {
  
  double dx = 0.;
  double dy = 0.;
  double dz = 0.;
  
  for (int i=0; i<mol->natoms; i++) {
    dx += mol->atoms[i].x * mol->atoms[i].tq;//indo;
    dy += mol->atoms[i].y * mol->atoms[i].tq;//indo;
    dz += mol->atoms[i].z * mol->atoms[i].tq;//indo;
 // cout<<mol->atoms[i].x<<" "<<mol->atoms[i].tq<<endl;
  }

  dx *= ang2au;
  dy *= ang2au;
  dz *= ang2au;

  mol->dipx = dx;
  mol->dipy = dy;
  mol->dipz = dz;
  
  cout<<"Transition dipole moment: "<<dx<<" x, "<<dy<<" y, "<<dz<<" z"<<endl;
  cout<<"Magnitude = "<<sqrt(dx*dx+dy*dy+dz*dz)/0.39345<<endl;
  return true;
}

/*** Calculate the transition charges for a given root ***/
bool calctqindo(Molecule *mol, int root) {
 double sum = 0.; 
  for (int atom=0; atom<mol->natoms; atom++) {

    mol->atoms[atom].tqindo = 0.;

    for (int b=0; b<mol->nbasis; b++) { //basis functions
      if (mol->nbasisatom[b] != atom+1) continue;
      for (int i=0; i<mol->nocc; i++) { //MO's ci's
        for (int j=mol->nocc; j<mol->nmo; j++) { //MO ci's
              mol->atoms[atom].tqindo += mol->ci[root*mol->nmo*mol->nmo+i+j*mol->nmo]
                  * mol->mos[i+b*mol->nmo] * mol->mos[j+b*mol->nmo];
        }//end unoccupied MO
      }//end occupied MO
    } //end NAO on atom of interest
    cout<<"atom "<<atom<<endl;
    mol->atoms[atom].tqindo *= sqrt(2);
    cout<<mol->nbasisatom[atom]<<" "<<mol->atoms[atom].type<<" "<<mol->atoms[atom].tqindo<<" "<<endl;
    sum += mol->atoms[atom].tqindo;
  } //end atoms
  
  cout<<"Net charge: "<<sum<<endl;
  return true;
}
/** Calculate WF overlap **/
bool calcOverlap(Molecule *mol, int root,Molecule *mol2, int root2) {
  double sum = 0.; 
  double s2,s3,s4;
  s2=0.;
  s3=0.;
  s4=0.;
  int active=30;
  
  //CTC Need to confirm this algorithm. Seems to work. 10-25-16
  for (int b=0; b<mol->nbasis; b++) { 
#pragma omp parallel for reduction (+:s2,s3,s4)
      for (int c=0; c<mol->nbasis; c++) { 
        for (int i=mol->nocc-active; i<mol->nocc; i++) {
          for (int j=mol->nocc; j< mol->nocc+active/*mol->nbasis*/; j++) {
            for (int k=mol->nocc-active; k<mol->nocc; k++) {
              for (int l=mol->nocc; l<mol->nocc+active /*mol->nbasis*/; l++) {
                if ((mol->cia[i+j*mol->nmo+root*mol->nmo*mol->nmo] == 0.)
                  && (mol2->cia[k+l*mol2->nmo+root2*mol2->nmo*mol2->nmo] == 0.))
                  continue;
              /*mol->atoms[atom].tq*///s2 += mol->ci[i+j*mol->nmo+root*mol->nmo*mol->nmo]
   /*                 * mol2->ci[k+l*mol2->nmo+root2*mol2->nmo*mol2->nmo]
                    * mol->mos[i+b*mol->nmo] 
                    * mol->mos[j+c*mol->nmo]
                    * mol2->mos[k+b*mol->nmo]
                    * mol2->mos[l+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];
              //if (mol->os) {*/
                  
                  
                  s3 += mol->cia[i+j*mol->nmo+root*mol->nmo*mol->nmo]
                    * mol2->cia[k+l*mol2->nmo+root2*mol2->nmo*mol2->nmo]
                    * mol->mos[i+b*mol->nmo] 
                    //* mol->mos[j+c*mol->nmo]
                    * mol2->mos[k+c*mol->nmo]
                    //* mol2->mos[l+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];

                /*mol->atoms[atom].tqb*/
                    s4 += mol->cib[i+j*mol->nmo+root*mol->nmo*mol->nmo]
                    * mol2->cib[k+l*mol2->nmo+root2*mol2->nmo*mol2->nmo]
                    * mol->mos[i+b*mol->nmo] 
                    * mol->mos[j+c*mol->nmo]
                    * mol2->mos[k+c*mol->nmo]
                    * mol2->mos[l+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];
          }//end unoccupied MO
        }//end occupied MO
      }//end all other AOs
    } //end NAO on atom of interest
  }
  cout<<b<<" out of "<<mol->nbasis<<" "<<s2<<" "<<s3<<" "<<s4<<" "<<s3+s4<<endl;
  }
    
  cout<<"Overlap sum= "<<s2<<" "<<s3<<" "<<s4<<endl;
  if (mol->os)
    mol->projection = s3+s4;
  else
    mol->projection = s2;
  return true;
}
/*** Calculate the transition charges for a given root ***/
bool calctq(Molecule *mol, int root) {
 double sum = 0.; 
  for (int atom=0; atom<mol->natoms; atom++) {
    double s2,s3,s4;
    s2=0.;
    s3=0.;
    s4=0.;
    mol->atoms[atom].tq = 0.;
    mol->atoms[atom].tqa = 0.;
    mol->atoms[atom].tqb = 0.;

    for (int b=0; b<mol->nbasis; b++) {
      if (mol->nbasisatom[b] != atom+1) continue;
      
#pragma omp parallel for reduction (+:s2,s3,s4)
      for (int c=0; c<mol->nbasis; c++) {
        //if (mol->nbasisatom[c] == atom+1) continue;
        
        for (int i=0; i<mol->nocc; i++) {
          for (int j=mol->nocc; j<mol->nbasis; j++) {
              /*mol->atoms[atom].tq*/s2 += mol->ci[i+j*mol->nmo+root*mol->nmo*mol->nmo]
                    * mol->mos[i+b*mol->nmo] * mol->mos[j+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];
              //if (mol->os) {
                /*mol->atoms[atom].tqa*/s3 += mol->cia[i+j*mol->nmo+root*mol->nmo*mol->nmo]
                    * mol->mos[i+b*mol->nmo] * mol->mos[j+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];
                /*mol->atoms[atom].tqb*/s4 += mol->cib[i+j*mol->nmo+root*mol->nmo*mol->nmo]
                    * mol->mos[i+b*mol->nmo] * mol->mos[j+c*mol->nmo]
                    * mol->overlapm[c+b*mol->nbasis];
              //}
          }//end unoccupied MO
        }//end occupied MO
      }//end all other AOs
    
    } //end NAO on atom of interest
    if (mol->os) {
      mol->atoms[atom].tqa = s3;
      mol->atoms[atom].tqb = s4;
      mol->atoms[atom].tq = s3+s4;
    } else
      mol->atoms[atom].tq = s2;

    /** Check on the sqrt(2), and if it's needed for open shell **/
    mol->atoms[atom].tq *= sqrt(2);
    mol->atoms[atom].tqa *= sqrt(2);
    mol->atoms[atom].tqb *= sqrt(2);
  cout<<atom<<" "<<mol->atoms[atom].type<<" "<<mol->atoms[atom].tq<<" "<<mol->atoms[atom].tqa<<" "<<mol->atoms[atom].tqb<<endl;
sum += mol->atoms[atom].tq;
  } //end atoms
  
  cout<<"Net charge: "<<sum<<endl;
  return true;
}

/*** Read in MO vectors from movec file ***/
/* str             - input, file name to read
 * mol            - input, Molecule object
 * mol->occupation - output, occupation numbers
 * mol->moeigenv   - output, MO eigenvalues
 * mol->mos        - output, MO in orthonormal basis
 */

bool readMOs(string str, Molecule *mol) {
  cout.precision(10);
  char tempc[1000];
  ifstream infile;
  infile.open(str.c_str());

  if (!infile.is_open()) {
    cout<<"Error opening MO file"<<endl;
    return false;
  } else {
    cout<<"Opening file: "<<str<<endl;
  }
  
  /* Skip right to the occupation numbers */
  getnlines(infile,tempc,14,500);
  mol->nocc = 0;
  mol->nuocc = 0;
  for (int i=0; i<mol->nmo; i++) {
    infile>>tempc;
    mol->occupation[i] = (int)atof(tempc);
    if (mol->occupation[i] == 0) {
      mol->nuocc++;
    } else {
      mol->nocc++;
    }
  }
cout<<"HOMO = "<<mol->nocc<<", # virt. orb. = "<<mol->nuocc<<endl;
  /* Read in eigenvalues */
  for (int i=0; i<mol->nmo; i++) {
    infile>>tempc;
    mol->moeigenv[i] = atof(tempc);
  }

  /* Read in MO vectors */
  for (int i=0; i<mol->nmo; i++)
    for (int j=0; j<mol->nbasis; j++) {
      infile>>tempc;
      mol->mos[i+j*mol->nmo] = atof(tempc);
    }
  return true;
}



/*** Function to parse an NWChem output file ***/
/* Returns true if successful, false otherwise */
/* 
 * str        - input, log file name
 * mol       - input, Molecule object
 *
 * mol->nroots    -output, number of tddft roots
 * mol->natoms    -output, number of atoms in molecule
 * mol->atom     -output, individual atom objects
 * mol->atom.type  - output, element
 * mol->atom.x,y,z - output, geometry
 * mol->atom.charge  output, charge on atom
 */

bool parseLog(string str, Molecule *mol) {
    
    int nci=0;
    ifstream infile;
    infile.precision(9);
    infile.open(str.c_str());
    cout<<"Opening file: "<<str.c_str()<<endl;

    char tempc[1000];
    int natoms=0;
    int nbasis=0;
    bool tddftstack = true;

    while(infile.getline(tempc, 1000)) {
      string temps(tempc);
      int current = infile.tellg();
      
/** Get Number of Roots in TDDFT **/
      if ((temps.compare(0,5,tddft_task_str1,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str2,0,5) == 0 ||
        temps.compare(0,5,tddft_task_str3,0,5) == 0) &&
        tddftstack) {
        cout<<"getting nroots"<<endl;
        while (temps.compare(0,3,"end",0,3) !=0 ) {
          if (temps.compare(0,6,"nroots")==0 ||
            temps.compare(0,6,"NROOTS")==0 ||
            temps.compare(0,6,"Nroots")==0 ||
            temps.compare(0,6,"NRoots")==0) {
          
            temps = strtok(tempc," ");
            temps = strtok(NULL," ");
            mol->setnroots(atoi(temps.c_str()));
            tddftstack = false;
          } //end nroots
          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        }
      /** create space for each root **/
      cout<<mol->nroots<<" nroots "<<endl;
      mol->roots = new Root[mol->nroots];
      for (int rt=0; rt<mol->nroots; rt++) {
        mol->roots[rt].ncia = 0;
        mol->roots[rt].ncib = 0;
      }

      } //end tddft


/** Geometry **/     
      if (temps.compare(0,10,geom_str,0,10) == 0) {
      cout<<"getting geometry"<<endl;
        /* Skip a few lines to get to the geometry */
        getnlines(infile,tempc,3,1000);

        /* Read geometry */
        //get position in file
/*        int pos = infile.tellg();
        while (temps != "Atomic") {
          infile>>ws;
          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          //get natoms
          if (temps != "Atomic") {
            natoms++;
          }
        }

        /* Make an molecule entry and declare natoms for the molecule */
//        mol->allocateMemAtoms(natoms);

        //return to top of geometry
//        infile.seekg(pos);
        infile.getline(tempc,1000);
        temps = strtok(tempc," ");
        int i=0;
        while(temps!= "Atomic") {
          //Assign atom numbers
          mol->atoms[i].num = atoi(temps.c_str());
          temps = strtok(NULL," ");

          //Get atom type
          mol->atoms[i].type = temps;

          //Get atom charge
          temps = strtok(NULL," ");
          mol->atoms[i].charge = atof(temps.c_str());

          //X-coordinate
          temps = strtok(NULL," ");
          mol->atoms[i].x = atof(temps.c_str());
          
          //Y-coordinate
          temps = strtok(NULL," ");
          mol->atoms[i].y = atof(temps.c_str());

          //Z-coordinate
          temps = strtok(NULL," ");
          mol->atoms[i].z = atof(temps.c_str());

          infile>>ws;
          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          i++;
        }
      } //end geometry

/** Number of AO basis functions **/
      
      if (temps.compare(0,40,nbasis_str,0,40)==0) {
        cout<<"getting nAO"<<endl;
        infile.getline(tempc,1000);
        temps = strtok(tempc,":");
        temps = strtok(NULL,":");

        //allocate memory for molecule's properties
        mol->allocateMem(atoi(temps.c_str()));
        for (int i=0; i<mol->nroots; i++) {
          mol->roots[i].ci = new double[mol->nmo];
        }

      } //end ao basis

/** Number linearly dependent vectors **/
    if (temps.compare(0,25,lindep_str,0,25)==0) {
      cout<<"getting nlindep vectors"<<endl;
      temps = strtok(tempc," ");
      for (int i=0; i<5; i++) temps = strtok(NULL," ");
      mol->allocateMemLindep(mol->nbasis,atoi(temps.c_str()));
    }


/** Basis Function Labels **/
      
      if (temps.compare(0,20,basislabel_str,0,20) == 0) {
      cout<<"getting basis function labels"<<endl;
        getnlines(infile,tempc,4,1000);
        temps = tempc;
        
        for (int i=0; i<mol->natoms; i++) 
          mol->atoms[i].nao = 0;

        while (temps.compare(0,10,overlap_str,1,10) != 0) {
          temps = tempc;
          temps = strtok(tempc," ");

          //get function number
          int num = atoi(temps.c_str());

          //get atom number
          temps = strtok(NULL," ");
          mol->nbasisatom[num-1] = atoi(temps.c_str());

          //increment number of AOs on atom
          mol->atoms[atoi(temps.c_str())-1].nao ++;
          
          //get atom element
          temps = strtok(NULL," ");
          mol->nbasisatomelements[num-1] = temps;

          //get type of AO
          temps = strtok(NULL," ");
          mol->nbasisatomorbitals[num-1] = temps;
          
          infile>>ws;
          infile.getline(tempc,1000);
          temps = tempc;
        } //end import basis functions

      } //end basis label

/** Basis function overlap matrix **/
      if (temps.compare(0,20,overlap_str,1,20) == 0) {
        cout<<"getting overlap matrix"<<endl;
        int nblocks = mol->nbasis/ncol;
        for (int block=0; block<nblocks; block++) {
          if (block == 0)
            getnlines(infile,tempc,4,1000);
          else
            getnlines(infile,tempc,3,1000);
          for (int i=0; i<mol->nbasis; i++) {
            temps = strtok(tempc," ");
            for (int j=0; j<ncol; j++) {
              temps = strtok(NULL," ");
              mol->overlapm[i+(j+ncol*block)*mol->nbasis]
                = atof(temps.c_str());
                //if (i==(j+ncol*block)) cout<<i<<" "<<mol->overlapm[i+(j+ncol*block)*mol->nbasis]<<endl;
            } // end col
              infile.getline(tempc,1000);
          }  // end rows
        } //end full blocks
        int remain = mol->nbasis%ncol;
        if (remain > 0) {
          getnlines(infile,tempc,3,1000);
        for (int i=0; i<mol->nbasis; i++) {
          temps = strtok(tempc," ");
          for (int j=0; j<remain; j++) {
            temps = strtok(NULL," ");
            mol->overlapm[i+(j+nblocks*ncol)*mol->nbasis] 
              = atof(temps.c_str());
            //if (i==189 && j==remain-1) cout<<mol->overlapm[i+(j+nblocks*ncol)*mol->nbasis]<<endl;
          } //end col
          infile.getline(tempc,1000);
        } //end row
        } //end remain
      } //end overlap matrix


/** TDDFT excited states **/
    if (temps.compare(0,10,tddft_str,0,10) == 0) {
      bool openshell = false;
      int dum = 1;
      cout<<"getting CI vectors"<<endl;
      mol->allocateMemTddft();
      getnlines(infile,tempc,3,1000);
      temps = tempc;
      if (temps.compare(0,5,"  <S2>",0,5) == 0) {
        getnlines(infile,tempc,2,1000);
        openshell = true;
        mol->os = true;
      } else
        infile.getline(tempc,1000);

      for (int root=0; root<mol->nroots; root++) {
        mol->roots[root].nci = 0;
        infile.getline(tempc,1000);
        temps = strtok(tempc," ");
        if (openshell) 
          dum = 3;
        else
          dum = 4;
        for (int i=0; i<dum; i++) temps = strtok(NULL," ");
        mol->excenergy[root] = atof(temps.c_str());
        if (openshell) {
          infile>>ws;
          infile.getline(tempc,1000);
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
          temps = strtok(NULL," ");
          mol->spin[root] = atof(temps.c_str());
        }
        if (openshell) 
          dum = 2;
        else
          dum = 2;

        getnlines(infile,tempc,dum,1000);
        
        for (int k=0; k<3; k++) {
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
          
          for (int j=0; j<3; j++) {
            for (int i=0; i<2; i++) temps = strtok(NULL," ");
              mol->transmoment[k+j*3+root*9] = atof(temps.c_str());
          }
          infile.getline(tempc,1000);
        }//end trans moment

        //get oscillator strength
        temps = strtok(tempc," ");
        for (int i=0; i<3; i++) temps=strtok(NULL," ");
        mol->oscstrength[root] = atof(temps.c_str());
        infile.getline(tempc,1000);
      
        infile.getline(tempc,1000);
        temps = tempc;
        
        //get CI coeffs
        //nrpa is for RPA
        /** Currently this is confirmed to work with 
         * closed shell tddft, and open shell cis
         **/
        int nrpa = 1;     
        double ycoeff;
        double sumci = 0.;
        double sumci2 = 0.;
        bool spina;

        
          nci=0;
        while (temps.compare(0,10,tddft_state_stop_str,0,10) != 0 &&
              temps.compare(0,10,"Target root",0,10) !=0 ) { //Start CI coeffs
          temps = strtok(tempc," ");
          temps = strtok(NULL," ");
          int row = atoi(temps.c_str())-1; //MO_from
          if (openshell) 
            dum = 5;
          else
            dum = 4;

          for (int i=0; i<dum; i++) temps = strtok(NULL," .");
          int col = atoi(temps.c_str())-1; //MO_to
          if (openshell) 
            dum = 1;
          else
            dum = 2;
          for (int i=0; i<dum; i++) temps = strtok(NULL," ");
          
          if (openshell) {
            if (temps.compare(0,4,"alpha",0,4)==0) {
              spina = true; //alpha
              mol->roots[root].ncia ++;
            } else {
              spina = false; //beta
              mol->roots[root].ncib ++;
            }

            for (int i=0; i<2; i++) temps = strtok(NULL," ");
          }
          mol->ci[row + col*mol->nmo + root*mol->nmo*mol->nmo] += atof(temps.c_str());
          mol->roots[root].ci[nci] = atof(temps.c_str());
            /** Alpha/Beta CI coefficients for unrestricted **/
            /** Alpha first, then beta, everything goes into roots.ci array of length ncia+ncib **/
            if (openshell && spina) {
              mol->cia[row + col*mol->nmo + root*mol->nmo*mol->nmo] += atof(temps.c_str());     
              mol->roots[root].ci[mol->roots[root].ncia-1] = atof(temps.c_str());
            } else if (openshell && (!spina)) {
              mol->cib[row + col*mol->nmo + root*mol->nmo*mol->nmo] += atof(temps.c_str());
              mol->roots[root].ci[
                  mol->roots[root].ncia
                  +mol->roots[root].ncib-1] = atof(temps.c_str());
            }

            sumci += atof(temps.c_str())*atof(temps.c_str());
          //}
          //CTC check this and make compatible with CIS 9-12-14
          /*if (nrpa%2 == 1 ) {
            mol->ci[row + col*mol->nmo + root*mol->nmo*mol->nmo] += atof(temps.c_str());
            sumci2 += atof(temps.c_str())*atof(temps.c_str());
          }
         */ nrpa++;
          nci++;
          infile>>ws;
          infile.getline(tempc,1000);
          temps=tempc;
        } //end CI coeffs for a given root
        mol->roots[root].nci = nci;
      } //end nroots
    }//end tddft
  } //end reading logfile
  
    return true;
}

bool getOverlapMatrix(string str, Molecule *mol) {
    
    /** gets AO overlap matrix from different geometries **/
    int nbas = 0;//mol->nbasis * 2;
    int nci=0;
    ifstream infile;
    infile.precision(9);
    infile.open(str.c_str());
    cout<<"Opening file: "<<str.c_str()<<endl;

    char tempc[1000];
    int natoms=0;
    int nbasis=0;
    bool tddftstack = true;

    while(infile.getline(tempc, 1000)) {
      string temps(tempc);
      /** Number of AO basis functions **/
      if (temps.compare(0,40,nbasis_str,0,40)==0) {
        cout<<"getting nAO for composite geometries"<<endl;
        infile.getline(tempc,1000);
        temps = strtok(tempc,":");
        temps = strtok(NULL,":");

        nbas = atoi(temps.c_str());
      } //end ao basis

      /** Basis function overlap matrix **/
      if (temps.compare(0,20,overlap_str,0,20) == 0) {
        cout<<"getting mutual overlap matrix"<<endl;
        int nblocks = nbas/ncol;
        for (int block=0; block<nblocks; block++) {
          if (block == 0)
            getnlines(infile,tempc,4,1000);
          else
            getnlines(infile,tempc,3,1000);
          for (int i=0; i<nbas; i++) {
            temps = strtok(tempc," ");
            for (int j=0; j<ncol; j++) {
              temps = strtok(NULL," ");
              mol->overlapm2[i+(j+ncol*block)*nbas]
                = atof(temps.c_str());
             // if (i==(j+ncol*block))
              //  cout<<i<<" "<<atof(temps.c_str())<<endl;
                //if (i==189 && j==5) cout<<mol->overlapm[i+(j+ncol*block)*mol->nbasis]<<endl;
            } // end col
              infile.getline(tempc,1000);
          }  // end rows
        } //end full blocks
        int remain = nbas%ncol;
        if (remain > 0) {
          getnlines(infile,tempc,3,1000);
        for (int i=0; i<nbas; i++) {
          temps = strtok(tempc," ");
          for (int j=0; j<remain; j++) {
            temps = strtok(NULL," ");
            mol->overlapm2[i+(j+nblocks*ncol)*nbas] 
              = atof(temps.c_str());
            //if (i==189 && j==remain-1) cout<<mol->overlapm[i+(j+nblocks*ncol)*mol->nbasis]<<endl;
          } //end col
          infile.getline(tempc,1000);
        } //end row
        } //end remain
      } //end overlap matrix
  }
//exit(0);
  /** Put mixed overlap matrix into overlapm **/
  /** NB: This is hard-coded for pmitet **/
  /** PMI-PMI overlap **/
  /** Actually this is not needed.
   * If this is done, the acceptor AO overlaps are
   * essentially 0 and the WF overlap is very small.
   * We want the character of each state, so we'll
   * take one of the geometry AO overlap matrices and 
   * use that
   */
//  for (int i=0; i<nbas; i++) {
    //mol->overlapm2[i+i*nbas] *= 1./(mol->overlapm2[i+i*nbas]);  
    //for (int j=0; j<nbas; j++) {
      //cout<<i<<" "<<j<<" "<<mol->overlapm2[i+j*nbas]<<endl;
    //}
//  }
/*  cout<<"Copying overlap matrices"<<endl;
  cout<<"composite nbas = "<<nbas<<endl;
  int nao_1 = 416;
  for (int i=0; i<416; i++) {
    for (int j=0; j<416; j++) {
      mol->overlapm[i+j*mol->nbasis] = mol->overlapm2[i+j*nbas];
    }
  }
  
  /** PMI-A1'A2' **/
/*  for (int i=0; i<416; i++) {
    for (int j=mol->nbasis; j<nbas; j++) {
      int col = j - 588;
      mol->overlapm[i+col*mol->nbasis] = mol->overlapm2[i+j*nbas];
    }
  }
  for (int i=416; i<mol->nbasis; i++) {
    for (int j=0; j<416; j++) {
      mol->overlapm[i+j*mol->nbasis] = mol->overlapm2[i+j*nbas];
    }
  }
  for (int i=416; i<mol->nbasis; i++) {
    for (int j=mol->nbasis; j<nbas; j++) {
      int col = j-588;
      mol->overlapm[i+col*mol->nbasis] = mol->overlapm2[i+j*nbas];
    }
  }
cout<<"Done with composite overlaps"<<endl;
  //for (int i=0; i<mol->nbasis; i++)
  //cout<<i<<" "<<mol->overlapm[i+i*mol->nbasis]<<endl;
*/
}

/** Print some stuff **/
void printStuff(Molecule *mol,string r1, string s1, string r2, string s2) {
  string filename = r1;
  filename.append("-");
  filename.append(s1);
  filename.append("-");
  filename.append(s2);
  filename.append(".dat");
  ofstream outfile;
  outfile.precision(16);
  outfile.open(filename.c_str(),ofstream::out | ofstream::app);
  
  /** print state and projection **/
  cout<<"Root "<<r2<<", projection= "<<mol->projection<<endl;
  outfile<<r2<<" "<<mol->projection<<endl;
}


#endif // UTIL_H

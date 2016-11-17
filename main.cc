/********************************************
 * Program to calculate atomic transition
 * densities using an AO expansion.
 *
 * CTC 7-31-2014
 *
 * For further information, see:
 * JPCB, 117, 2032, 2013
 * JPCC, 114, 20834, 2010
 * *******************************************/

#include "main.h"
#include <vector>
using namespace std;

int main(int argc, char **argv) {
  
  //open input file
  comfile.open(argv[1]);


  //parse input file
  char tempc[1000];
  string temps;

  /* 
   * Get calculation type
   * "transden" - get ao transition densities
   * "overlap"  - get overlap of two chromophore
   *              transition densities
   */
  
    /* Allocate memory for molecule object */
    Molecule molecule, molecule2;

    /* get natoms */
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL," :");

    /* Make a molecule entry and declare natoms for the molecules */
    molecule.allocateMemAtoms(atoi(temps.c_str()));
    molecule2.allocateMemAtoms(atoi(temps.c_str()));
    
    /* get nwchem log file name for molecule 1*/
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL," :");
    string temps3 = temps;

    /* get MO vectors from movec file for molecule 1*/
    /* File must be generated using the
     * mov2asc utility provided with NWChem
     */
    /* get movecs file name */
    comfile.getline(tempc,1000);
    temps = strtok(tempc," :");
    temps=strtok(NULL," :");
    string temps2 = temps;

    /* get nwchem log file name for molecule 2*/
    comfile.getline(tempc,1000);
    temps = strtok(tempc,":");
    temps = strtok(NULL," :");
    string temps4 = temps;
cout<<"getting tddft from "<<temps2<<", and "<<temps4<<endl;
    /* get MO vectors from movec file for molecule 2*/
    /* File must be generated using the
     * mov2asc utility provided with NWChem
     */
    /* get movecs file name */
    comfile.getline(tempc,1000);
    temps = strtok(tempc," :");
    temps=strtok(NULL," :");
    string temps5 = temps;


    /* Excitation level of theory */
    comfile.getline(tempc,1000);
    temps = strtok(tempc, " :");
    temps = strtok(NULL, " ");
    molecule.excMethod = temps;
    cout<<"Excitation method: "<<temps<<endl;
     
    /* read in which root to consider for molecule 1*/
    comfile.getline(tempc,1000);
    temps = strtok(tempc,": ");
    temps = strtok(NULL,": ");
    string sroot1 = temps;
    int nroot = atoi(temps.c_str())-1;   
    
    /* read in which root to consider for molecule 2*/
    comfile.getline(tempc,1000);
    temps = strtok(tempc,": ");
    temps = strtok(NULL,": ");
    string sroot2 = temps;
    int nroot2 = atoi(temps.c_str())-1;   
    cout<<"Using root "<<nroot<<" for mol 1, and "<<nroot2<<" for mol 2."<<endl;
    
    /* read in chromophore separation for cluster 1 */
    comfile.getline(tempc,1000);
    temps = strtok(tempc,": ");
    temps = strtok(NULL,": ");
    string sep1 = temps;   
    molecule.separation = atof(temps.c_str());

    /* read in chromophore separation for cluster 1 */
    comfile.getline(tempc,1000);
    temps = strtok(tempc,": ");
    temps = strtok(NULL,": ");
    string sep2 = temps;   
    molecule2.separation = atof(temps.c_str());

    /* read in log file for AO overlaps */
    comfile.getline(tempc,1000);
    temps = strtok(tempc,": ");
    temps = strtok(NULL,": ");
    string aofile = temps;
    
    /* Parse Log file */
    parseLog(temps3,&molecule);
    parseLog(temps4,&molecule2);
    cout<<"Done parsing log files"<<endl;

    /** Get cross AO overlap **/
    //getOverlapMatrix(aofile,&molecule);

    /* read in mo's */
    readMOs(temps2,&molecule);
    readMOs(temps5,&molecule2);

    cout<<"Done reading MOs"<<endl;
   
   /** Calculate WF overlap **/
   calcOverlap2(&molecule,nroot,&molecule2,nroot2);

    /* Calculate transition dipole moment */
    //calcdip(&molecule);

    /* Print stuff */
    //printStuff(&molecule,sroot1,sep1,sroot2,sep2);

  return 0;
}

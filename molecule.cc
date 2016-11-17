#include "molecule.h"

Molecule::Molecule()
{
}

void Molecule::allocateMem(const int nb) {
  this->nbasis = nb;
  this->nbasisatom = new int[nb];
  this->nbasisatomelements = new string[nb];
  this->nbasisatomorbitals = new string[nb];
  this->overlapm = new double[nb*nb];
  this->overlapm2 = new double[4*nb*nb]; //cross-mol ao overlap
  this->occupation = new int[nb];
  this->moeigenv = new double[nb];
  this->mos = new double[nb*nb];
  this->nlindep = 0;
  this->nmo = nb;//= mol->nbasis;
  this->projection = 0.;
}

void Molecule::allocateMemLindep(const int nb, const int nlindep) {
  //delete[] this->mos;
  this->nlindep = nlindep;
  this->nmo = nb-nlindep;
  //this->mos = new double[this->nmo*nb];
}

void Molecule::allocateMemAtoms(const int na) {
  this->natoms = na;
  this->atoms = new Atom[na];
  this->activeCharges = new double[na];
}

void Molecule::allocateMemTddft() {
  this->os = false;
  int ntrans = 9*this->nroots;
  this->spin = new double[this->nroots];
  this->excenergy = new double[this->nroots];
  this->transmoment = new double[ntrans];
  this->oscstrength = new double[this->nroots];
  this->ci = new double[this->nroots*this->nmo*this->nmo];
  this->cia = new double[this->nroots*this->nmo*this->nmo];
  this->cib = new double[this->nroots*this->nmo*this->nmo];
  for (int i=0; i<this->nroots*this->nmo*this->nmo; i++) {
    this->ci[i] = 0.;
    this->cia[i] = 0.;
    this->cib[i] = 0.;
    }
}

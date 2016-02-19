/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include "writer.h"
#include "view.h"
#include <libmints/mints.h>
#include <psi4-dec.h>
#include <physconst.h>

#include <cstdio>
#include <utility>
#include <algorithm>
#include "libparallel/ParallelPrinter.h"
using namespace std;
using namespace psi;
using namespace boost;

GradientWriter::GradientWriter(boost::shared_ptr<Molecule> mol, const Matrix& grad)
    : molecule_(mol), gradient_(grad)
{
}

void GradientWriter::write(const std::string &filename)
{
   boost::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));
   int i;


    printer->Printf("%-59.59s %-10.10s%-9.9s\n",
            molecule_->name().c_str(),
            "(wfn)",
            "(dertype)");

    printer->Printf("%5d%20.10lf\n", molecule_->natom(), Process::environment.globals["CURRENT ENERGY"]);

    for (i=0; i<molecule_->natom(); ++i) {
        printer->Printf("%20.10lf%20.10lf%20.10lf%20.10lf\n",
                double(molecule_->Z(i)), molecule_->x(i), molecule_->y(i), molecule_->z(i));
    }

    for (i=0; i<molecule_->natom(); ++i) {
        printer->Printf("                    %20.10lf%20.10lf%20.10lf\n",
                gradient_(i, 0), gradient_(i, 1), gradient_(i, 2));
    }
}

MoldenWriter::MoldenWriter(boost::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{

}
void MoldenWriter::writeNO(const std::string &filename, boost::shared_ptr<Matrix> Na, boost::shared_ptr<Matrix> Nb, boost::shared_ptr<Vector> Oa, boost::shared_ptr<Vector> Ob)
{
    //Same as MO Writer below 
    boost::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));

    int atom;

    printer->Printf("[Molden Format]\n");
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    // Print the molecule to molden
    printer->Printf("[Atoms] (AU)\n");
    for (atom=0; atom<mol.natom(); ++atom) {
        Vector3 coord = mol.xyz(atom);
        printer->Printf("%-2s  %2d  %3d   %20.12f %20.12f %20.12f\n",
                mol.symbol(atom).c_str(), atom+1, static_cast<int>(mol.Z(atom)), coord[0], coord[1], coord[2]);
    }

    // Dump the basis set using code adapted from psi2molden
    printer->Printf("[GTO]\n");

    // For each atom
    for (atom=0; atom<mol.natom(); ++atom) {
        printer->Printf("  %d 0\n", atom+1);

        // Go through all the shells on this center
        for (int shell=0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            const GaussianShell& gs = basisset.shell(overall_shell);

            printer->Printf(" %c%5d  1.00\n", gs.amchar(), gs.nprimitive());

            for (int prim=0; prim<gs.nprimitive(); ++prim) {
                printer->Printf("%20.10f %20.10f\n", gs.exp(prim), gs.original_coef(prim));
            }
        }

        // An empty line separates atoms
        printer->Printf("\n");
    }
    /* Natural Orbital Transformation to AO basis 
     *N  (mo x no)
     *1st Half Transform
     *N' (so x no) = C (so x mo) x N (mo x no)
     *Fully transformed 
     *N'' (ao x no) = S (ao x no) x N'(so x no)
     */
    //setup 
    // get the "S" transformation matrix, ao by so
    boost::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));
    SharedMatrix aotoso = pl->aotoso();
    //get C's
    SharedMatrix Ca = wavefunction_->Ca();
    SharedMatrix Cb = wavefunction_->Cb();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension nmo = Ca->colspi();
    const Dimension nos = Na->colspi();
    //New N's
    SharedMatrix Naprime(new Matrix("Na' ", sos, nos));
    SharedMatrix Nbprime(new Matrix("Nb' ", sos, nos));
    // do N' = C x N 
    Naprime->gemm(false, false, 1.0, Ca,Na, 0.0);
    Nbprime->gemm(false, false, 1.0, Cb,Na, 0.0);
    //Fully transformed
    SharedMatrix NaFT(new Matrix("NaFT", aos,nos));
    SharedMatrix NbFT(new Matrix("NbFT", aos,nos));
    // do N'' = S x N'
    NaFT->gemm(false,false,1.0,aotoso,Naprime,0.0);
    NbFT->gemm(false,false,1.0,aotoso,Nbprime,0.0);
    
    // The order Molden expects
    //     P: x, y, z
    //    5D: D 0, D+1, D-1, D+2, D-2
    //    6D: xx, yy, zz, xy, xz, yz
    //
    //    7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
    //   10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    //
    //    9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
    //   15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
    //        xxyy xxzz yyzz xxyz yyxz zzxy
    // Since Molden doesn't handle higher than g we'll just leave them be.
    int molden_cartesian_order[][15] = {
        { 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // p
        { 0, 3, 4, 1, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // d
        { 0, 4, 5, 3, 9, 6, 1, 8, 7, 2, 0, 0, 0, 0, 0 },    // f
        { 0, 3, 4, 9, 12, 10, 5, 13, 14, 7, 1, 6, 11, 8, 2} // g
    };

    int nirrep = NaFT->nirrep();
    Dimension countpi(nirrep);
    Dimension zeropi(nirrep);
    Dimension ncartpi(nirrep);

    for(int i = 0; i < basisset.nshell(); i++) {
        int am = basisset.shell(i).am();

        int ncart = basisset.shell(i).nfunction();
        if((am == 1 && basisset.has_puream()) || (am > 1 && am < 5 && basisset.shell(i).is_cartesian())) {
            for (int h=0; h<nirrep; ++h)
                ncartpi[h] = ncart;

            View block_a(NaFT, ncartpi, NaFT->colspi(), countpi, zeropi);
            View block_b(NbFT, ncartpi, NbFT->colspi(), countpi, zeropi);

            SharedMatrix temp_a = block_a();
            SharedMatrix temp_b = block_b();

            for( int j =0; j < ncart; j++) {
                for (int h=0; h < NaFT->nirrep(); ++h) {
                    for (int k=0; k<NaFT->coldim(h); ++k) {
                        // outfile->Printf( "am %d\n, from %d to %d\n", am, j, countpi[h] + molden_cartesian_order[am-1][j]);
                        NaFT->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_a->get(h, j, k));
                        NaFT->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_b->get(h, j, k));
                    }
                }
            }
        }

        for (int h=0; h<nirrep; ++h)
            countpi[h] += ncart;
    }

    if (basisset.has_puream()) {
        // Tell Molden to use spherical.  5d implies 5d and 7f.
        printer->Printf("[5D]\n[9G]\n\n");
    }
    CharacterTable ct = mol.point_group()->char_table();

    // Dump MO's to the molden file
    printer->Printf("[MO]\n");

    std::vector<std::pair<double, std::pair<int, int> > > mos;

    // do alpha's
    bool SameOcc = true;
    for (int h=0; h<wavefunction_->nirrep(); ++h) {
        for (int n=0; n<wavefunction_->nmopi()[h]; ++n) {
            mos.push_back(make_pair(Oa->get(h, n), make_pair(h, n)));
            if(fabs(Oa->get(h,n) - Ob->get(h,n)) > 1e-10)
                SameOcc = false;
        }
    }
    std::sort(mos.begin(), mos.end());

    for (int i=0; i<(int)mos.size(); ++i) {
        int h = mos[i].second.first;
        int n = mos[i].second.second;

        printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
        printer->Printf(" Ene= %20.10f\n", Oa->get(h, n));
        printer->Printf(" Spin= Alpha\n");
        if(Na == Nb && Oa == Ob && SameOcc)
            printer->Printf(" Occup= %7.4lf\n", Oa->get(h,n)+Ob->get(h,n));
        else
            printer->Printf(" Occup= %7.4lf\n", Oa->get(h,n));
        for (int so=0; so<wavefunction_->nso(); ++so)
            printer->Printf("%3d %20.12lf\n", so+1, NaFT->get(h, so, n));
    }

    // do beta's
    mos.clear();
    if (Na != Nb || Oa != Ob || !SameOcc) {
        for (int h=0; h<wavefunction_->nirrep(); ++h) {
            for (int n=0; n<wavefunction_->nmopi()[h]; ++n) {
                mos.push_back(make_pair(Ob->get(h, n), make_pair(h, n)));
            }
        }
        std::sort(mos.begin(), mos.end());

        for (int i=0; i<(int)mos.size(); ++i) {
            int h = mos[i].second.first;
            int n = mos[i].second.second;

            printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
            printer->Printf(" Ene= %20.10lf\n", Ob->get(h, n));
            printer->Printf(" Spin= Beta\n");
            printer->Printf(" Occup= %7.4lf\n", Ob->get(h,n));
            for (int so=0; so<wavefunction_->nso(); ++so)
                printer->Printf("%3d %20.12lf\n", so+1, NbFT->get(h, so, n));
        }
    }


}


void MoldenWriter::write(const std::string &filename, boost::shared_ptr<Matrix> Ca, boost::shared_ptr<Matrix> Cb, boost::shared_ptr<Vector> Ea, boost::shared_ptr<Vector> Eb, boost::shared_ptr<Vector> OccA, boost::shared_ptr<Vector> OccB)
{
    boost::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));

    int atom;

    printer->Printf("[Molden Format]\n");

    // Get the molecule for ease
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    //    basisset.print_detail();

    // Print the molecule to molden
    printer->Printf("[Atoms] (AU)\n");
    for (atom=0; atom<mol.natom(); ++atom) {
        Vector3 coord = mol.xyz(atom);
        printer->Printf("%-2s  %2d  %3d   %20.12f %20.12f %20.12f\n",
                mol.symbol(atom).c_str(), atom+1, static_cast<int>(mol.Z(atom)), coord[0], coord[1], coord[2]);
    }

    // Dump the basis set using code adapted from psi2molden
    printer->Printf("[GTO]\n");

    // For each atom
    for (atom=0; atom<mol.natom(); ++atom) {
        printer->Printf("  %d 0\n", atom+1);

        // Go through all the shells on this center
        for (int shell=0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            const GaussianShell& gs = basisset.shell(overall_shell);

            printer->Printf(" %c%5d  1.00\n", gs.amchar(), gs.nprimitive());

            for (int prim=0; prim<gs.nprimitive(); ++prim) {
                printer->Printf("%20.10f %20.10f\n", gs.exp(prim), gs.original_coef(prim));
            }
        }

        // An empty line separates atoms
        printer->Printf("\n");
    }

    // Convert Ca & Cb
    boost::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));
    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension nmo = Ca->colspi();

    SharedMatrix Ca_ao_mo(new Matrix("Ca AO x MO", aos, nmo));
    SharedMatrix Cb_ao_mo(new Matrix("Cb AO x MO", aos, nmo));

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    //    aotoso->print();
    //    Ca_ao_mo->print();
    //    Cb_ao_mo->print();

    // The order Molden expects
    //     P: x, y, z
    //    5D: D 0, D+1, D-1, D+2, D-2
    //    6D: xx, yy, zz, xy, xz, yz
    //
    //    7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
    //   10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    //
    //    9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
    //   15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
    //        xxyy xxzz yyzz xxyz yyxz zzxy
    // Since Molden doesn't handle higher than g we'll just leave them be.
    int molden_cartesian_order[][15] = {
        { 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // p
        { 0, 3, 4, 1, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // d
        { 0, 4, 5, 3, 9, 6, 1, 8, 7, 2, 0, 0, 0, 0, 0 },    // f
        { 0, 3, 4, 9, 12, 10, 5, 13, 14, 7, 1, 6, 11, 8, 2} // g
    };

    int nirrep = Ca_ao_mo->nirrep();
    Dimension countpi(nirrep);
    Dimension zeropi(nirrep);
    Dimension ncartpi(nirrep);

    for(int i = 0; i < basisset.nshell(); i++) {
        int am = basisset.shell(i).am();

        int ncart = basisset.shell(i).nfunction();
        if((am == 1 && basisset.has_puream()) || (am > 1 && am < 5 && basisset.shell(i).is_cartesian())) {
            for (int h=0; h<nirrep; ++h)
                ncartpi[h] = ncart;

            View block_a(Ca_ao_mo, ncartpi, Ca_ao_mo->colspi(), countpi, zeropi);
            View block_b(Cb_ao_mo, ncartpi, Cb_ao_mo->colspi(), countpi, zeropi);

            SharedMatrix temp_a = block_a();
            SharedMatrix temp_b = block_b();

            for( int j =0; j < ncart; j++) {
                for (int h=0; h < Ca_ao_mo->nirrep(); ++h) {
                    for (int k=0; k<Ca_ao_mo->coldim(h); ++k) {
                        // outfile->Printf( "am %d\n, from %d to %d\n", am, j, countpi[h] + molden_cartesian_order[am-1][j]);
                        Ca_ao_mo->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_a->get(h, j, k));
                        Cb_ao_mo->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_b->get(h, j, k));
                    }
                }
            }
        }

        for (int h=0; h<nirrep; ++h)
            countpi[h] += ncart;
    }

    if (basisset.has_puream()) {
        // Tell Molden to use spherical.  5d implies 5d and 7f.
        printer->Printf("[5D]\n[9G]\n\n");
    }
    CharacterTable ct = mol.point_group()->char_table();

    // Dump MO's to the molden file
    printer->Printf("[MO]\n");

    std::vector<std::pair<double, std::pair<int, int> > > mos;

    // do alpha's
    bool SameOcc = true;
    for (int h=0; h<wavefunction_->nirrep(); ++h) {
        for (int n=0; n<wavefunction_->nmopi()[h]; ++n) {
            mos.push_back(make_pair(Ea->get(h, n), make_pair(h, n)));
            if(fabs(OccA->get(h,n) - OccB->get(h,n)) > 1e-10)
                SameOcc = false;
        }
    }
    std::sort(mos.begin(), mos.end());

    for (int i=0; i<(int)mos.size(); ++i) {
        int h = mos[i].second.first;
        int n = mos[i].second.second;

        printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
        printer->Printf(" Ene= %20.10f\n", Ea->get(h, n));
        printer->Printf(" Spin= Alpha\n");
        if(Ca == Cb && Ea == Eb && SameOcc)
            printer->Printf(" Occup= %7.4lf\n", OccA->get(h,n)+OccB->get(h,n));
        else
            printer->Printf(" Occup= %7.4lf\n", OccA->get(h,n));
        for (int so=0; so<wavefunction_->nso(); ++so)
            printer->Printf("%3d %20.12lf\n", so+1, Ca_ao_mo->get(h, so, n));
    }

    // do beta's
    mos.clear();
    if (Ca != Cb || Ea != Eb || !SameOcc) {
        for (int h=0; h<wavefunction_->nirrep(); ++h) {
            for (int n=0; n<wavefunction_->nmopi()[h]; ++n) {
                mos.push_back(make_pair(Eb->get(h, n), make_pair(h, n)));
            }
        }
        std::sort(mos.begin(), mos.end());

        for (int i=0; i<(int)mos.size(); ++i) {
            int h = mos[i].second.first;
            int n = mos[i].second.second;

            printer->Printf(" Sym= %s\n", ct.gamma(h).symbol());
            printer->Printf(" Ene= %20.10lf\n", Eb->get(h, n));
            printer->Printf(" Spin= Beta\n");
            printer->Printf(" Occup= %7.4lf\n", OccB->get(h,n));
            for (int so=0; so<wavefunction_->nso(); ++so)
                printer->Printf("%3d %20.12lf\n", so+1, Cb_ao_mo->get(h, so, n));
        }
    }


}

NBOWriter::NBOWriter(boost::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{


}

void NBOWriter::write(const std::string &filename)
{
    int pure_order[][7] = {
        { 1, 0, 0, 0, 0, 0, 0},      // s
        { 101, 102, 103, 0, 0, 0, 0}, // p
        // z2  xz   yz  x2-y2 xy
        { 255, 252, 253, 254, 251, 0, 0}, // d
        //z(z2-r2), x(z2-r2), y(z2-r2) z(x2-y2), xyz, x(x2-y2), y(x2-y2)
        { 351, 352, 353, 354, 355, 356, 357 } //f
    };

    MintsHelper helper(wavefunction_->basisset(), wavefunction_->options(), 0);
    SharedMatrix sotoao = helper.petite_list()->sotoao();
    boost::shared_ptr<OutFile> printer(new OutFile(filename,APPEND));


    //Get the basis set and molecule from the wavefuntion
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    //NBO can only handle up to f functions
    if( basisset.max_am () > 3)
    {
        throw PSIEXCEPTION("NBO cannot handle angular momentum above f functions. \n");
    }
    //print $GENNBO section of file
    //BOHR indicates atomic units for the coordinates; now ANG but not sure about keyword
    //OPEN indicates that we'll provide separate alpha and beta matrices
    printer->Printf(" $GENNBO NATOMS = %d NBAS = %d BODM ", mol.natom(), basisset.nbf());

    //To make this more user-friendly in the case of RHF wavefunctions...
    bool open_shell = (wavefunction_->nalpha() != wavefunction_->nbeta());
    if(open_shell)
        printer->Printf(" OPEN $END\n");
    else
        printer->Printf(" $END\n");

    //print NBO section of file47; user can modify this to suit their needs
    printer->Printf(" $NBO       $END\n");

    //Now print out the molecule
    printer->Printf(" $COORD\n");
    printer->Printf(" GENNBO expects one comment line here. So, here's a comment line.\n");
    for( int i =0; i< mol.natom(); i++)
    {
        //the second mol.Z() should be modified when pseudopotentials are implemented
        printer->Printf( "%2d  %2d  %20.12f %20.12f %20.12f\n",
                static_cast<int>(mol.Z(i)), static_cast<int>(mol.Z(i)),
                mol.x(i)*pc_bohr2angstroms, mol.y(i)*pc_bohr2angstroms,
                mol.z(i)*pc_bohr2angstroms);
    }
    printer->Printf( " $END\n");


    //To form the BASIS and CONTRACT sections, we need some memory
    int nshells = basisset.nshell(); //Total number of shells
    int nprim = basisset.nprimitive(); //Total number of primitives
    Vector centers(basisset.nbf());
    Vector labels(basisset.nbf());
    Vector components(nshells); //Functions per shell
    Vector angmom(nshells); //Angular momentum of shell
    Vector nprimitives(nshells); //Primitives per shell
    Vector exponents(nprim); //Exponents of primitives
    //Coefficient matrix; first row is S, second P, third D, fourth F
    Matrix coefficient(4, nprim);
    coefficient.zero();
    int fnindex = 0;
    int primindex = 0;

    //Loop over all the shells
    for( int i =0; i < nshells; i++)
    {
        const GaussianShell& gshell = basisset.shell(i);
        int nfns = gshell.nfunction(); //get number of functions in shell
        components.set(0, i, nfns);
        int angm = gshell.am(); //get angular momentum of shell
        angmom.set(0, i, angm);
        for( int j = 0; j< nfns; j++)
        {
            centers.set (0, fnindex, gshell.ncenter());
            if(gshell.is_pure()) {
                //outfile->Printf( "fnindex %d pure_order[%d][%d] %d\n", fnindex, angm, j, pure_order[angm][j]);
                labels.set (0, fnindex, pure_order[angm][j]);
            }
            else
                labels.set (0, fnindex, angm*100+j+1);
            fnindex++;
        }
        int nshellprim = gshell.nprimitive();
        nprimitives.set (0, i, nshellprim);
        for( int k =0; k < nshellprim; k++)
        {
            exponents.set(0, primindex, gshell.exp(k));
            coefficient.set (0, angm, primindex, gshell.coef(k));
            primindex++;
        }
    }

    // now, we print out the basis section
    printer->Printf(" $BASIS\n");
    // CENTER section
    for(int i=0; i<basisset.nbf(); i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n  CENTER =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)centers.get(0, i)+1);
    }

    //The LABEL section
    for( int i =0; i < basisset.nbf(); i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n   LABEL =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)labels.get(0, i));
    }
    printer->Printf( "\n $END\n");

    //The CONTRACT heading
    printer->Printf( " $CONTRACT\n");
    printer->Printf( "  NSHELL = %6d\n", nshells);
    printer->Printf( "    NEXP = %6d\n", nprim);

    // List of the number of functions per shell
    for(int i=0; i < nshells; i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n   NCOMP =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)components.get(0, i));
    }
    // List the number of primitives per shell
    for(int i=0; i < nshells; i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n   NPRIM =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", (int)nprimitives.get(0, i));
    }
    // location of the first exponent for each shell
    int ptr = 1;
    for( int i=0; i < nshells; i++)
    {
        if((i+1)%10 == 1) {
            if(i==0)
                printer->Printf("\n    NPTR =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( " %6d", ptr);
        ptr += nprimitives.get(0, i);
    }
    // exponents
    for( int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n     EXP =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", exponents.get(0, i));
    }
    // coefficients for s functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CS =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 0, i));
    }
    // coefficients for p functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CP =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 1, i));
    }
    // coefficients for d functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CD =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 2, i));
    }
    // coefficients for f functions
    for(int i=0; i<nprim; i++)
    {
        if((i+1)%4 == 1) {
            if(i==0)
                printer->Printf("\n      CF =");
            else
                printer->Printf("\n          ");
        }
        printer->Printf( "%15.6E", coefficient.get (0, 3, i));
    }
    printer->Printf( "\n $END");

    //Matrix transformation information we'll need
    int nbf = basisset.nbf ();

    //Now we need the overlap matrix in the AO basis
    SharedMatrix overlap = helper.ao_overlap();
    //Print overlap matrix
    printer->Printf( "\n $OVERLAP");
    for(int i=0; i<nbf; i++)
    {
        for(int j=0; j<nbf; j++)
        {
            if(((nbf*i+j+1)%5) == 1)
                printer->Printf("\n  ");
            printer->Printf( "%15.6E", overlap->get (0, i, j));
        }
    }
    printer->Printf( "\n $END");

    //Alpha Density Matrix
    SharedMatrix soalphadens = wavefunction_->Da();
    SharedMatrix alphadens(new Matrix(nbf, nbf));
    alphadens->remove_symmetry (soalphadens, sotoao);
    //Beta density
    SharedMatrix betadens(new Matrix(nbf, nbf));
    SharedMatrix sobetadens = wavefunction_->Db();
    betadens->remove_symmetry (sobetadens, sotoao);
    //Now print the density matrix
    printer->Printf( "\n $DENSITY");
    if(wavefunction_->same_a_b_dens ())
    {
        SharedMatrix density(new Matrix(nbf, nbf));
        density->copy (alphadens);
        density->add (betadens);
        for( int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                if(((nbf*i+j+1)%5)==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", density->get(0, i, j));
            }
        }
    }
    else
    {
        int count = 0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                count++;
                if(count%5 == 1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphadens->get (0, i, j));
            }
        }
        for(int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                count++;
                if(count%5 ==0)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", betadens->get (0, i, j));
            }
        }
    }
    printer->Printf( "\n $END");


    // alpha Fock matrix
    SharedMatrix alphasofock = wavefunction_->Fa();
    SharedMatrix alphafock(new Matrix(nbf, nbf));
    alphafock->remove_symmetry (alphasofock, sotoao);
    // print the Fock matrix
    printer->Printf( "\n $FOCK");
    if(wavefunction_->same_a_b_dens ())
    {
        for(int i = 0; i < nbf; i++)
        {
            for(int j = 0; j < nbf; j++)
            {
                if(((nbf*i+j+1)%5)==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphafock->get (0, i, j));
            }
        }
    }

    else
    {
        // beta Fock
        SharedMatrix betafock(new Matrix(nbf, nbf));
        SharedMatrix betasofock = wavefunction_->Fb();
        betafock->remove_symmetry(betasofock, sotoao);
        int count=0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                count++;
                if(count%5 == 1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphafock->get (0, i, j));
            }
        }
        for(int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                count++;
                if(count%5 == 1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", betafock->get (0, i, j));
            }
        }
    }
    printer->Printf( "\n $END");

    //Alpha AO->MO transformation
    SharedMatrix soalphac = wavefunction_->Ca();
    const Dimension aos = helper.petite_list()->AO_basisdim();
    const Dimension nmo = wavefunction_->Ca()->colspi();
    SharedMatrix alphac(new Matrix("Ca AO x MO", aos, nmo));
    alphac->gemm(true, false, 1.00, sotoao, soalphac, 0.00);

    printer->Printf( "\n $LCAOMO");
    if(wavefunction_->same_a_b_orbs ())
    {
        for(int i = 0; i < nbf; i++)
        {
            for(int j = 0; j < nbf; j++)
            {
                if(((nbf*i+j+1)%5)==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphac->get(0, i, j));
            }
        }
    }

    else
    {
        //Beta AO->MO transformation
        SharedMatrix betac(new Matrix(nbf, nbf));
        SharedMatrix sobetac = wavefunction_->Cb();
        betac->gemm(true, false, 1.00, sotoao, sobetac, 0.00);

        //Print the AO->MO coefficients
        int count = 0;
        for(int i=0; i<nbf; i++)
        {
            for(int j=0; j<nbf; j++)
            {
                count++;
                if(count%5==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E", alphac->get (0, i, j));
            }
        }
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                count++;
                if(count%5 ==1)
                    printer->Printf( "\n  ");
                printer->Printf( "%15.6E ", betac->get (0,   i, j));
            }
        }
    }
    printer->Printf( "\n $END\n");

}


MOWriter::MOWriter(boost::shared_ptr<Wavefunction> wavefunction,Options&options)
    : wavefunction_(wavefunction), options_(options)
{
}

void MOWriter::write()
{
    // what kind of reference?
    bool isrestricted = true;
    if (options_.get_str("REFERENCE") == "UHF")  isrestricted = false;
    if (options_.get_str("REFERENCE") == "UKS")  isrestricted = false;
    if (options_.get_str("REFERENCE") == "CUHF") isrestricted = false;

    // Get the molecule for ease
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule & mol = *basisset.molecule().get();

    // Convert Ca & Cb
    // make copies
    Matrix Ca(wavefunction_->Ca());
    Matrix Cb(wavefunction_->Cb());
    Vector& Ea = *wavefunction_->epsilon_a().get();
    Vector& Eb = *wavefunction_->epsilon_b().get();

    boost::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));

    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();
    const Dimension mos = wavefunction_->nmopi();

    SharedMatrix Ca_ao_mo(new Matrix("Ca AO x MO", aos, mos));
    SharedMatrix Cb_ao_mo(new Matrix("Cb AO x MO", aos, mos));

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    int nirrep = Ca_ao_mo->nirrep();

    // order orbitals in terms of energy:
    int minorb;
    nmo = mos.sum();

    map  = new int[nmo];
    bool * skip = new bool[nmo];
    for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
    for (int orb = 0; orb < nmo; orb++) {

        int count = 0;
        double minen = 1.0e9;
        for (int h = 0; h < nirrep; h++) {
            for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {

                if ( skip[count] ) {
                    count++;
                    continue;
                }
                if ( Ea.get(h,n) <= minen ) {
                    minen = Ea.get(h,n);
                    minorb = count;
                }

                count++;

            }
        }
        map[ orb ] = minorb;
        skip[ minorb ] = true;
    }

    // reorder orbitals:
    nso = wavefunction_->nso();
    eps = new double[nmo];
    sym = new int[nmo];
    occ = new int[nmo];
    Ca_pointer = new double[nmo * nso];
    for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;

    int count = 0;
    int extra = isrestricted ? 1 : 0;
    for (int h = 0; h < nirrep; h++) {
        double ** Ca_old = Ca_ao_mo->pointer(h);
        for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
            occ[ count ] = n < ( wavefunction_->doccpi()[h] + wavefunction_->soccpi()[h] ) ? 1 : 0;
            occ[ count ] += n <  wavefunction_->doccpi()[h] ? extra : 0;
            eps[ count ] = Ea.get(h,n);
            sym[ count ] = h;
            for (int mu = 0; mu < nso; mu++) {
                Ca_pointer[mu*nmo + count] = Ca_old[mu][n];
            }
            count++;

        }
    }

    // dump to output file
    outfile->Printf("\n");
    if ( isrestricted )
        outfile->Printf("  ==> Molecular Orbitals <==\n");
    else
        outfile->Printf("  ==> Alpha-Spin Molecular Orbitals <==\n");
    outfile->Printf("\n");

    write_mos(mol);

    // now for beta spin
    if ( !isrestricted ) {


        // order orbitals in terms of energy
        for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
        for (int orb = 0; orb < nmo; orb++) {

            int count = 0;
            double minen = 1.0e9;
            for (int h = 0; h < nirrep; h++) {
                for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {

                    if ( skip[count] ) {
                        count++;
                        continue;
                    }

                    if ( Eb.get(h,n) <= minen ) {
                        minen = Eb.get(h,n);
                        minorb = count;
                    }
                    count++;

                }
            }
            map[ orb ] = minorb;
            skip[ minorb ] = true;
        }


        // reorder orbitals:
        for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;
        count = 0;
        for (int h = 0; h < nirrep; h++) {
            double ** Ca_old = Cb_ao_mo->pointer(h);
            for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
                occ[ count ] = n < wavefunction_->doccpi()[h] ? 1 : 0;
                eps[ count ] = Eb.get(h,n);
                sym[ count ] = h;
                for (int mu = 0; mu < nso; mu++) {
                    Ca_pointer[mu*nmo + count] = Ca_old[mu][n];
                }
                count++;

            }
        }

        // dump to output file
        outfile->Printf("\n");
        outfile->Printf("  ==> Beta-Spin Molecular Orbitals <==\n");
        outfile->Printf("\n");
        write_mos(mol);
    }

    delete skip;
    delete occ;
    delete sym;
    delete eps;
    delete Ca_pointer;
}

void MOWriter::write_mos(Molecule & mol){

    CharacterTable ct = mol.point_group()->char_table();

    // print mos (5 columns)
    int ncols = 5;
    int ncolsleft = nmo % ncols;
    int nrows = (nmo - ncolsleft ) / ncols;

    // print the full rows:
    int count = 0;
    for (int i = 0; i < nrows; i++) {

        // print blank space
        outfile->Printf("     ");
        // print mo number
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13d",count+j+1);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao number
            outfile->Printf("%5i",mu+1);
            for (int j = 0; j < ncols; j++){

                outfile->Printf("%13.7lf",Ca_pointer[ mu*nmo + map[count + j] ]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
        // print energy
        outfile->Printf(" Ene ");
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13.7lf",eps[ map[count + j] ]);
        }
        outfile->Printf("\n");
        // print symmetry
        outfile->Printf(" Sym ");
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13s",ct.gamma(sym[map[count+j]]).symbol());
        }
        outfile->Printf("\n");
        // print occupancy
        outfile->Printf(" Occ ");
        for (int j = 0; j < ncols; j++){
            outfile->Printf("%13d",occ[map[count+j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        outfile->Printf("\n");
        count+=ncols;
    }

    // print the partial rows:
    if ( ncolsleft > 0 ) {

        // print blank space
        outfile->Printf("     ");
        // print mo number
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13d",count+j+1);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao number
            outfile->Printf("%5i",mu+1);
            for (int j = 0; j < ncolsleft; j++){
                outfile->Printf("%13.7lf",Ca_pointer[ mu*nmo + map[count + j] ]);
            }
            outfile->Printf("\n");
        }
        outfile->Printf("\n");
        // print energy
        outfile->Printf(" Ene ");
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13.7lf",eps[ map[count + j] ]);
        }
        outfile->Printf("\n");
        outfile->Printf(" Sym ");
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13s",ct.gamma(sym[map[count+j]]).symbol());
        }
        outfile->Printf("\n");
        // print occupancy
        outfile->Printf(" Occ ");
        for (int j = 0; j < ncolsleft; j++){
            outfile->Printf("%13d",occ[map[count+j]]);
        }
        outfile->Printf("\n");
        outfile->Printf("\n");

    }
}


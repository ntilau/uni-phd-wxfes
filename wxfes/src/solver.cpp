#include "solver.h"

#include <smumps_c.h>
#include <dmumps_c.h>
#include <cmumps_c.h>
#include <zmumps_c.h>

#include <cfloat>

static double shape_grad[] = { -1.0, 1.0, 0.0, 0.0,
                               -1.0, 0.0, 1.0, 0.0,
                               -1.0, 0.0, 0.0, 1.0
                               };
// EDGE
static std::vector<double> wq1_felippa = { 1.184634425280945e-01,
                                           2.393143352496832e-01,
                                           2.844444444444444e-01,
                                           2.393143352496832e-01,
                                           1.184634425280945e-01
                                         };
static std::vector<double> xq1_felippa = { 4.691007703066802e-02,
                                           2.307653449471585e-01,
                                           5.000000000000000e-01,
                                           7.692346550528415e-01,
                                           9.530899229693319e-01
                                         };
static std::vector<double> xq1_kythe = { 3.084018393622490e-03,
                                         1.975436564598987e-02,
                                         5.577038356387148e-02,
                                         1.127016653792583e-01,
                                         1.894485266313868e-01,
                                         2.828781253265987e-01,
                                         3.883066567855166e-01,
                                         5.000000000000000e-01,
                                         6.116933432144834e-01,
                                         7.171218746734013e-01,
                                         8.105514733686132e-01,
                                         8.872983346207417e-01,
                                         9.442296164361286e-01,
                                         9.802456343540101e-01,
                                         9.969159816063775e-01
                                       };
static std::vector<double> wq1_kythe = { 8.500859814970131e-03,
                                         2.580164149853987e-02,
                                         4.646359765756227e-02,
                                         6.720762762189211e-02,
                                         8.575595456819569e-02,
                                         1.003142646884945e-01,
                                         1.095784292007938e-01,
                                         1.127552498991034e-01,
                                         1.095784292007938e-01,
                                         1.003142646884945e-01,
                                         8.575595456819569e-02,
                                         6.720762762189211e-02,
                                         4.646359765756227e-02,
                                         2.580164149853987e-02,
                                         8.500859814970131e-03
                                       };

// TRIANGLE
static std::vector<double> wq2_dunavant1 = { 5.000000000000000e-001 };
static std::vector<std::vector<double> > xq2_dunavant1 = {
    { 3.333333333333333e-001, 3.333333333333333e-001 }
};
static std::vector<double> wq2_strang1 = { 1.666666666666667e-001,
                                           1.666666666666667e-001,
                                           1.666666666666667e-001
                                         };
static std::vector<std::vector<double> > xq2_strang1 = {
    { 6.666666666666666e-001, 1.666666666666667e-001 },
    { 1.666666666666667e-001, 6.666666666666666e-001 },
    { 1.666666666666667e-001, 1.666666666666667e-001 }
};
static std::vector<double> wq2_strang5 = { 5.497587182766100e-002,
                                           5.497587182766100e-002,
                                           5.497587182766100e-002,
                                           1.116907948390055e-001,
                                           1.116907948390055e-001,
                                           1.116907948390055e-001
                                         };
static std::vector<std::vector<double> > xq2_strang5 = {
    { 8.168475729804590e-001, 9.157621350977101e-002 },
    { 9.157621350977101e-002, 8.168475729804590e-001 },
    { 9.157621350977101e-002, 9.157621350977101e-002 },
    { 1.081030181680700e-001, 4.459484909159650e-001 },
    { 4.459484909159650e-001, 1.081030181680700e-001 },
    { 4.459484909159650e-001, 4.459484909159650e-001 }
};
static std::vector<double> wq2_strang8 = { 1.029752523804435e-001,
                                           1.029752523804435e-001,
                                           1.029752523804435e-001,
                                           3.184570714311150e-002,
                                           3.184570714311150e-002,
                                           3.184570714311150e-002,
                                           3.184570714311150e-002,
                                           3.184570714311150e-002,
                                           3.184570714311150e-002
                                         };

static std::vector<std::vector<double> > xq2_strang8 = {
    { 1.249495032332320e-001, 4.375252483833840e-001 },
    { 4.375252483833840e-001, 1.249495032332320e-001 },
    { 4.375252483833840e-001, 4.375252483833840e-001 },
    { 7.971126518600710e-001, 1.654099273898410e-001 },
    { 7.971126518600710e-001, 3.747742075008800e-002 },
    { 1.654099273898410e-001, 7.971126518600710e-001 },
    { 1.654099273898410e-001, 3.747742075008800e-002 },
    { 3.747742075008800e-002, 7.971126518600710e-001 },
    { 3.747742075008800e-002, 1.654099273898410e-001 }
};
static std::vector<double> wq2_strang9 = { 2.542245318510350e-002,
                                           2.542245318510350e-002,
                                           2.542245318510350e-002,
                                           5.839313786318950e-002,
                                           5.839313786318950e-002,
                                           5.839313786318950e-002,
                                           4.142553780918700e-002,
                                           4.142553780918700e-002,
                                           4.142553780918700e-002,
                                           4.142553780918700e-002,
                                           4.142553780918700e-002,
                                           4.142553780918700e-002
                                         };
static std::vector<std::vector<double> > xq2_strang9 = {
    { 8.738219710169960e-001, 6.308901449150201e-002 },
    { 6.308901449150201e-002, 8.738219710169960e-001 },
    { 6.308901449150201e-002, 6.308901449150201e-002 },
    { 2.492867451709100e-001, 5.014265096581790e-001 },
    { 5.014265096581790e-001, 2.492867451709100e-001 },
    { 2.492867451709100e-001, 2.492867451709100e-001 },
    { 6.365024991213990e-001, 3.103524510337850e-001 },
    { 6.365024991213990e-001, 5.314504984481600e-002 },
    { 3.103524510337850e-001, 6.365024991213990e-001 },
    { 3.103524510337850e-001, 5.314504984481600e-002 },
    { 5.314504984481600e-002, 6.365024991213990e-001 },
    { 5.314504984481600e-002, 3.103524510337850e-001 }
};

// TETRAHEDRON
static std::vector<double> wq3_keast1 = { 4.166666666666666e-002,
                                          4.166666666666666e-002,
                                          4.166666666666666e-002,
                                          4.166666666666666e-002
                                        };
static std::vector<std::vector<double> > xq3_keast1 = {
    { 5.854101966249685e-001, 1.381966011250105e-001, 1.381966011250105e-001 },
    { 1.381966011250105e-001, 1.381966011250105e-001, 1.381966011250105e-001 },
    { 1.381966011250105e-001, 1.381966011250105e-001, 5.854101966249685e-001 },
    { 1.381966011250105e-001, 5.854101966249685e-001, 1.381966011250105e-001 }
};
static std::vector<double> wq3_keast4 = { -1.315555555555555e-002,
                                          7.622222222222217e-003,
                                          7.622222222222217e-003,
                                          7.622222222222217e-003,
                                          7.622222222222217e-003,
                                          2.488888888888888e-002,
                                          2.488888888888888e-002,
                                          2.488888888888888e-002,
                                          2.488888888888888e-002,
                                          2.488888888888888e-002,
                                          2.488888888888888e-002
                                        };
static std::vector<std::vector<double> > xq3_keast4 = {
    { 2.500000000000000e-001, 2.500000000000000e-001, 2.500000000000000e-001 },
    { 7.857142857142857e-001, 7.142857142857140e-002, 7.142857142857140e-002 },
    { 7.142857142857140e-002, 7.142857142857140e-002, 7.142857142857140e-002 },
    { 7.142857142857140e-002, 7.142857142857140e-002, 7.857142857142857e-001 },
    { 7.142857142857140e-002, 7.857142857142857e-001, 7.142857142857140e-002 },
    { 1.005964238332008e-001, 3.994035761667992e-001, 3.994035761667992e-001 },
    { 3.994035761667992e-001, 1.005964238332008e-001, 3.994035761667992e-001 },
    { 3.994035761667992e-001, 3.994035761667992e-001, 1.005964238332008e-001 },
    { 3.994035761667992e-001, 1.005964238332008e-001, 1.005964238332008e-001 },
    { 1.005964238332008e-001, 3.994035761667992e-001, 1.005964238332008e-001 },
    { 1.005964238332008e-001, 1.005964238332008e-001, 3.994035761667992e-001 }
};
static std::vector<double> wq3_keast7 = { 6.653791709694645e-003,
                                          6.653791709694645e-003,
                                          6.653791709694645e-003,
                                          6.653791709694645e-003,
                                          1.679535175886776e-003,
                                          1.679535175886776e-003,
                                          1.679535175886776e-003,
                                          1.679535175886776e-003,
                                          9.226196923942399e-003,
                                          9.226196923942399e-003,
                                          9.226196923942399e-003,
                                          9.226196923942399e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003,
                                          8.035714285714283e-003
                                        };
static std::vector<std::vector<double> > xq3_keast7 = {
    { 3.561913862225449e-001, 2.146028712591517e-001, 2.146028712591517e-001 },
    { 2.146028712591517e-001, 2.146028712591517e-001, 2.146028712591517e-001 },
    { 2.146028712591517e-001, 2.146028712591517e-001, 3.561913862225449e-001 },
    { 2.146028712591517e-001, 3.561913862225449e-001, 2.146028712591517e-001 },
    { 8.779781243961660e-001, 4.067395853461134e-002, 4.067395853461134e-002 },
    { 4.067395853461134e-002, 4.067395853461134e-002, 4.067395853461134e-002 },
    { 4.067395853461134e-002, 4.067395853461134e-002, 8.779781243961660e-001 },
    { 4.067395853461134e-002, 8.779781243961660e-001, 4.067395853461134e-002 },
    { 3.298632957317306e-002, 3.223378901422757e-001, 3.223378901422757e-001 },
    { 3.223378901422757e-001, 3.223378901422757e-001, 3.223378901422757e-001 },
    { 3.223378901422757e-001, 3.223378901422757e-001, 3.298632957317306e-002 },
    { 3.223378901422757e-001, 3.298632957317306e-002, 3.223378901422757e-001 },
    { 2.696723314583159e-001, 6.366100187501753e-002, 6.366100187501753e-002 },
    { 6.366100187501753e-002, 2.696723314583159e-001, 6.366100187501753e-002 },
    { 6.366100187501753e-002, 6.366100187501753e-002, 2.696723314583159e-001 },
    { 6.030056647916491e-001, 6.366100187501753e-002, 6.366100187501753e-002 },
    { 6.366100187501753e-002, 6.030056647916491e-001, 6.366100187501753e-002 },
    { 6.366100187501753e-002, 6.366100187501753e-002, 6.030056647916491e-001 },
    { 6.366100187501753e-002, 2.696723314583159e-001, 6.030056647916491e-001 },
    { 2.696723314583159e-001, 6.030056647916491e-001, 6.366100187501753e-002 },
    { 6.030056647916491e-001, 6.366100187501753e-002, 2.696723314583159e-001 },
    { 6.366100187501753e-002, 6.030056647916491e-001, 2.696723314583159e-001 },
    { 2.696723314583159e-001, 6.366100187501753e-002, 6.030056647916491e-001 },
    { 6.030056647916491e-001, 2.696723314583159e-001, 6.366100187501753e-002 }
};

solver::solver(mdl_core& _mdl) :
    mdl(_mdl) {
    compute_dof_num();
    compute_label_map();
    if (strcmp(mdl.frm.type.data(), "EM_E_FD") == 0)
        analyze_em_e_fd();
    else if (strcmp(mdl.frm.type.data(), "E_V_STAT") == 0)
        analyze_e_v_stat();
    else if (strcmp(mdl.frm.type.data(), "H_A_STAT") == 0)
        analyze_h_a_stat();
    else
        std::cout << "Formulation not yet implemented\n";
}

solver::~solver() {
    //dtor
}

void solver::solve_dmumps(sp_mat<double>& matrix, sp_mat<double>& rhs) {
    std::vector<std::vector<double> >().swap(mdl.frm.sol_real);
    //mapping();
    size_t nrows = matrix.size();
    size_t a_nnz = matrix.nnz();
    std::cout << "N = " << nrows << ", NNZ = " << a_nnz << " ";
    DMUMPS_STRUC_C* idz = new DMUMPS_STRUC_C;
    idz->job = -1; // init
    idz->par = 1;
    idz->sym = (MUMPS_INT) 0;
    idz->comm_fortran = -987654; // use_comm_world
    dmumps_c(idz);
    idz->n = (MUMPS_INT) nrows;
    idz->nz = (MUMPS_INT) a_nnz;
    idz->irn = new MUMPS_INT[idz->nz];
    idz->jcn = new MUMPS_INT[idz->nz];
    idz->a = new DMUMPS_REAL[idz->nz];
    size_t idx = 0;
    for (typename sp_mat<double>::iterator itr = matrix.begin();
            itr != matrix.end(); itr++) {
        for (typename sp_vec<double>::iterator itc = itr->begin();
                itc != itr->end(); itc++) {
            idz->a[idx] = (DMUMPS_REAL) itc->val;
            idz->irn[idx] = (MUMPS_INT) (itr - matrix.begin()) + 1;
            idz->jcn[idx++] = (MUMPS_INT) itc->idx + 1;
        }
    }
    matrix.clear();
    idz->icntl[0] = -1; //-1 6
    idz->icntl[1] = -1; //-1 6
    idz->icntl[2] = -1; //-1 1
    idz->icntl[3] = 0; //0 6
    idz->icntl[19] = 1; // sparse right hand side
    idz->nrhs = (MUMPS_INT) rhs.size();
    idz->lrhs = (MUMPS_INT) nrows;
    idz->rhs = new DMUMPS_REAL[idz->lrhs * idz->nrhs];
    idz->nz_rhs = (MUMPS_INT) rhs.nnz();
    idz->rhs_sparse = new DMUMPS_REAL[idz->nz_rhs];
    idz->irhs_sparse = new MUMPS_INT[idz->nz_rhs];
    idz->irhs_ptr = new MUMPS_INT[idz->nrhs + 1];
    idx = 0;
    size_t i = 0;
    for (typename sp_mat<double>::iterator itr = rhs.begin(); itr != rhs.end();
            itr++) {
        idz->irhs_ptr[i++] = (MUMPS_INT) idx + 1;
        for (typename sp_vec<double>::iterator itc = itr->begin();
                itc != itr->end(); itc++) {
            if (itc->idx < nrows) {
                idz->irhs_sparse[idx] = (MUMPS_INT) itc->idx + 1;
                idz->rhs_sparse[idx++] = (DMUMPS_REAL) itc->val;
            }
        }
    }
    idz->irhs_ptr[idz->nrhs] = (MUMPS_INT) (idz->nz_rhs + 1);
    idz->job = 6; // factorization and solution
    dmumps_c(idz);
    std::cout << "-- " << idz->info[14] << " MB\n";
    //raw_sol.resize(raw_to_cmpct.size(), std::vector < T > (idz->nrhs, T(0)));
    mdl.frm.sol_real.resize(idz->nrhs); // , std::vector<double>(idz->lrhs)
    size_t shift;
    for (size_t col = 0; col < idz->nrhs; col++) {
        shift = col * idz->lrhs;
        for (size_t row = 0; row < idz->lrhs; row++) {
            mdl.frm.sol_real[col].push_back(idz->rhs[shift + row]);
        }
    }
    idz->job = -2; // end
    dmumps_c(idz);
    delete idz->a;
    delete idz->rhs;
    delete idz->irn;
    delete idz->jcn;
    delete idz->rhs_sparse;
    delete idz->irhs_sparse;
    delete idz->irhs_ptr;
    delete idz;
}

void solver::solve_zmumps(sp_mat<std::complex<double> >& matrix,
                          sp_mat<std::complex<double> >& rhs) {
    std::vector<std::vector<std::complex<double>> >().swap(mdl.frm.sol_cmplx);
    //mapping();
    size_t nrows = matrix.size();
    size_t a_nnz = matrix.nnz();
    std::cout << "N = " << nrows << ", NNZ = " << a_nnz << " ";
    ZMUMPS_STRUC_C* idz = new ZMUMPS_STRUC_C;
    idz->job = -1; // init
    idz->par = 1;
    idz->sym = (MUMPS_INT) 0;
    idz->comm_fortran = -987654; // use_comm_world
    zmumps_c(idz);
    idz->n = (MUMPS_INT) nrows;
    idz->nz = (MUMPS_INT) a_nnz;
    idz->irn = new MUMPS_INT[idz->nz];
    idz->jcn = new MUMPS_INT[idz->nz];
    idz->a = new ZMUMPS_COMPLEX[idz->nz];
    size_t idx = 0;
    for (typename sp_mat<std::complex<double> >::iterator itr = matrix.begin();
            itr != matrix.end(); itr++) {
        for (typename sp_vec<std::complex<double> >::iterator itc =
                    itr->begin(); itc != itr->end(); itc++) {
            idz->a[idx].r = (ZMUMPS_REAL) itc->val.real();
            idz->a[idx].i = (ZMUMPS_REAL) itc->val.imag();
            idz->irn[idx] = (MUMPS_INT) (itr - matrix.begin()) + 1;
            idz->jcn[idx++] = (MUMPS_INT) itc->idx + 1;
        }
    }
    matrix.clear();
    idz->icntl[0] = -1; //-1 6
    idz->icntl[1] = -1; //-1 6
    idz->icntl[2] = -1; //-1 1
    idz->icntl[3] = 0; //0 6
    idz->icntl[19] = 1; // sparse right hand side
    idz->nrhs = (MUMPS_INT) rhs.size();
    idz->lrhs = (MUMPS_INT) nrows;
    idz->rhs = new ZMUMPS_COMPLEX[idz->lrhs * idz->nrhs];
    idz->nz_rhs = (MUMPS_INT) rhs.nnz();
    idz->rhs_sparse = new ZMUMPS_COMPLEX[idz->nz_rhs];
    idz->irhs_sparse = new MUMPS_INT[idz->nz_rhs];
    idz->irhs_ptr = new MUMPS_INT[idz->nrhs + 1];
    idx = 0;
    size_t i = 0;
    for (typename sp_mat<std::complex<double> >::iterator itr = rhs.begin();
            itr != rhs.end(); itr++) {
        idz->irhs_ptr[i++] = (MUMPS_INT) idx + 1;
        for (typename sp_vec<std::complex<double> >::iterator itc =
                    itr->begin(); itc != itr->end(); itc++) {
            if (itc->idx < nrows) {
                idz->irhs_sparse[idx] = (MUMPS_INT) itc->idx + 1;
                idz->rhs_sparse[idx].r = (ZMUMPS_REAL) itc->val.real();
                idz->rhs_sparse[idx++].i = (ZMUMPS_REAL) itc->val.imag();
            }
        }
    }
    idz->irhs_ptr[idz->nrhs] = (MUMPS_INT) (idz->nz_rhs + 1);
    idz->job = 6; // factorization and solution
    zmumps_c(idz);
    std::cout << "-- " << idz->info[14] << " MB\n";
    //raw_sol.resize(raw_to_cmpct.size(), std::vector < T > (idz->nrhs, T(0)));
    mdl.frm.sol_cmplx.resize(idz->nrhs); // , std::vector<double>(idz->lrhs)
    size_t shift;
    for (size_t col = 0; col < idz->nrhs; col++) {
        shift = col * idz->lrhs;
        for (size_t row = 0; row < idz->lrhs; row++) {
            mdl.frm.sol_cmplx[col].push_back(
                std::complex<double>(idz->rhs[shift + row].r,
                                     idz->rhs[shift + row].i));
        }
    }
    idz->job = -2; // end
    zmumps_c(idz);
    delete idz->a;
    delete idz->rhs;
    delete idz->irn;
    delete idz->jcn;
    delete idz->rhs_sparse;
    delete idz->irhs_sparse;
    delete idz->irhs_ptr;
    delete idz;
}

void solver::compute_label_map() {
    std::cout << "Mapping labels of materials";
    for (unsigned int i = 0; i < mdl.frm.mtrls.size(); i++) {
        mtrl_map[mdl.frm.mtrls[i].label] = i;
        mdl.frm.mtrls[i].tetras.clear();
        mdl.frm.mtrls[i].faces.clear();
    }
    std::cout << ", boundaries";
    for (unsigned int i = 0; i < mdl.frm.bcs.size(); i++) {
        bc_map[mdl.frm.bcs[i].label] = i;
        mdl.frm.bcs[i].faces.clear();
        mdl.frm.bcs[i].edges.clear();
    }
    if (mdl.msh.n_tetras == 0) {
        std::cout << ", in 2D model\n";
        for (size_t i = 0; i < mdl.msh.n_faces; i++) {
            mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[i]]].faces.push_back(i);
        }
        for (size_t i = 0; i < mdl.msh.n_edges; i++) {
            if (mdl.msh.edg_lab[i] != 0)
                mdl.frm.bcs[bc_map[mdl.msh.edg_lab[i]]].edges.push_back(i);
        }
    } else {
        std::cout << ", in 3D model\n";
        for (size_t i = 0; i < mdl.msh.n_tetras; i++) {
            mdl.frm.mtrls[mtrl_map[mdl.msh.tet_lab[i]]].tetras.push_back(i);
        }
        for (size_t i = 0; i < mdl.msh.n_faces; i++) {
            if (mdl.msh.fac_lab[i] != 0)
                mdl.frm.bcs[bc_map[mdl.msh.fac_lab[i]]].faces.push_back(i);
        }
    }
}

void solver::compute_dof_num() {
    switch (mdl.frm.p) {
    case 1:
        dof_num.hgrad = mdl.msh.n_nodes;
        dof_num.hcurl = mdl.msh.n_edges;
        dof_num.hdiv = 0;
        break;
    case 2:
        dof_num.hgrad = mdl.msh.n_nodes + mdl.msh.n_edges;
        dof_num.hcurl = 2 * (mdl.msh.n_edges + mdl.msh.n_faces);
        dof_num.hdiv = 0;
        break;
    case 3:
        dof_num.hgrad = mdl.msh.n_nodes + 2 * mdl.msh.n_edges + mdl.msh.n_faces;
        dof_num.hcurl = 3 * (mdl.msh.n_edges + mdl.msh.n_tetras)
                        + 6 * mdl.msh.n_faces;
        dof_num.hdiv = 0;
        break;
    }
    std::cout << "Raw DoFs number = \n"
              << "Hgrad = " << dof_num.hgrad << "\n"
              << "Hcurl = " << dof_num.hcurl << "\n"
              << "Hdiv = " << dof_num.hdiv << "\n";
}

std::vector<size_t> solver::get_dofs(type s, type e, size_t id) {
    std::vector<size_t> dofs;
    switch (s) {
    case HGRAD:
        switch (e) {
        case EDGE:
            dofs = mdl.msh.edg_nodes[id];
            if (mdl.frm.p > 1)
                dofs.push_back(mdl.msh.n_nodes + id);
            if (mdl.frm.p > 2)
                dofs.push_back(mdl.msh.n_nodes + mdl.msh.n_edges + id);
            break;
        case TRIA:
            dofs = mdl.msh.fac_nodes[id];
            if (mdl.frm.p > 1)
                for (unsigned int i = 0; i < 3; i++)
                    dofs.push_back(mdl.msh.n_nodes + mdl.msh.fac_edges[id][i]);
            if (mdl.frm.p > 2) {
                for (unsigned int i = 0; i < 3; i++)
                    dofs.push_back(
                        mdl.msh.n_nodes + mdl.msh.n_edges
                        + mdl.msh.fac_edges[id][i]);
                dofs.push_back(mdl.msh.n_nodes + 2 * mdl.msh.n_edges + id);
            }
            break;
        case TETRA:
            dofs = mdl.msh.tet_nodes[id];
            if (mdl.frm.p > 1)
                for (unsigned int i = 0; i < 6; i++)
                    dofs.push_back(mdl.msh.n_nodes + mdl.msh.tet_edges[id][i]);
            if (mdl.frm.p > 2) {
                for (unsigned int i = 0; i < 6; i++)
                    dofs.push_back(
                        mdl.msh.n_nodes + mdl.msh.n_edges
                        + mdl.msh.tet_edges[id][i]);
                for (unsigned int i = 0; i < 4; i++)
                    dofs.push_back(
                        mdl.msh.n_nodes + 2 * mdl.msh.n_edges
                        + mdl.msh.tet_faces[id][i]);
            }
            break;
        }
        break;
    case HCURL:
        switch (e) {
        case EDGE:
            dofs.push_back(id);
            if (mdl.frm.p > 1)
                dofs.push_back(mdl.msh.n_edges + id);
            if (mdl.frm.p > 2)
                dofs.push_back(2 * mdl.msh.n_edges + 2 * mdl.msh.n_faces + id);
            break;
        case TRIA:
            dofs = mdl.msh.fac_edges[id];
            if (mdl.frm.p > 1) {
                for (unsigned int i = 0; i < 3; i++)
                    dofs.push_back(mdl.msh.n_edges + mdl.msh.fac_edges[id][i]);
                dofs.push_back(2 * mdl.msh.n_edges + id);
                dofs.push_back(2 * mdl.msh.n_edges + mdl.msh.n_faces + id);
            }
            if (mdl.frm.p > 2) {
                for (unsigned int i = 0; i < 3; i++)
                    dofs.push_back(
                        2 * mdl.msh.n_edges + 2 * mdl.msh.n_faces
                        + mdl.msh.fac_edges[id][i]);
                dofs.push_back(3 * mdl.msh.n_edges + 2 * mdl.msh.n_faces + id);
                dofs.push_back(3 * mdl.msh.n_edges + 3 * mdl.msh.n_faces + id);
                dofs.push_back(3 * mdl.msh.n_edges + 4 * mdl.msh.n_faces + id);
                dofs.push_back(3 * mdl.msh.n_edges + 5 * mdl.msh.n_faces + id);
            }
            break;
        case TETRA:
            dofs = mdl.msh.tet_edges[id];
            if (mdl.frm.p > 1) {
                for (unsigned int i = 0; i < 6; i++)
                    dofs.push_back(mdl.msh.n_edges + mdl.msh.tet_edges[id][i]);
                for (unsigned int i = 0; i < 4; i++) {
                    dofs.push_back(
                        2 * mdl.msh.n_edges + mdl.msh.tet_faces[id][i]);
                    dofs.push_back(
                        2 * mdl.msh.n_edges + mdl.msh.n_faces
                        + mdl.msh.tet_faces[id][i]);
                }
            }
            if (mdl.frm.p > 2) {
                for (unsigned int i = 0; i < 6; i++)
                    dofs.push_back(
                        2 * mdl.msh.n_edges + 2 * mdl.msh.n_faces
                        + mdl.msh.tet_edges[id][i]);
                for (unsigned int i = 0; i < 4; i++) {
                    dofs.push_back(
                        3 * mdl.msh.n_edges + 2 * mdl.msh.n_faces
                        + mdl.msh.tet_faces[id][i]);
                    dofs.push_back(
                        3 * mdl.msh.n_edges + 3 * mdl.msh.n_faces
                        + mdl.msh.tet_faces[id][i]);
                    dofs.push_back(
                        3 * mdl.msh.n_edges + 4 * mdl.msh.n_faces
                        + mdl.msh.tet_faces[id][i]);
                    dofs.push_back(
                        3 * mdl.msh.n_edges + 5 * mdl.msh.n_faces
                        + mdl.msh.tet_faces[id][i]);
                }
                dofs.push_back(3 * mdl.msh.n_edges + 6 * mdl.msh.n_faces + id);
                dofs.push_back(
                    3 * mdl.msh.n_edges + 6 * mdl.msh.n_faces
                    + mdl.msh.n_tetras + id);
                dofs.push_back(
                    3 * mdl.msh.n_edges + 6 * mdl.msh.n_faces
                    + 2 * mdl.msh.n_tetras + id);
            }
            break;
        }
        break;
    case HDIV:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            break;
        case TETRA:
            break;
        }
        break;
    }
    return dofs;
}

std::pair<double, std::vector<std::vector<double> > > solver::get_jac(
    std::vector<std::vector<double> > geo) {
    double det;
    std::vector<std::vector<double> > jac, inv;
    unsigned int n = geo.size();
    jac.resize(n - 1, std::vector<double>(n - 1, 0.0));
    for (unsigned int i = 0; i < n - 1; i++)
        for (unsigned int j = 0; j < n - 1; j++)
            for (unsigned int k = 0; k < n; k++)
                jac[i][j] += shape_grad[4 * i + k] * geo[k][j];
    inv.resize(n - 1, std::vector<double>(n - 1, 0.0));
    switch (n - 1) {
    case 1:
        det = std::sqrt(
                  std::pow(std::abs(geo[1][0] - geo[0][0]), 2)
                  + std::pow(std::abs(geo[1][1] - geo[0][1]), 2)
                  + std::pow(std::abs(geo[1][2] - geo[0][2]), 2));
        inv[0][0] = 1.0;
        break;
    case 2:
        det = std::abs(+jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]);
        inv[0][0] = jac[1][1] / det;
        inv[1][0] = -jac[0][1] / det;
        inv[0][1] = -jac[1][0] / det;
        inv[1][1] = jac[0][0] / det;
        break;
    case 3:
        det =
            std::abs(+ jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2])
                     - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
                     + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]));
        inv[0][0] = (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) / det;
        inv[1][0] = -(jac[0][1] * jac[2][2] - jac[0][2] * jac[2][1]) / det;
        inv[2][0] = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) / det;
        inv[0][1] = -(jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) / det;
        inv[1][1] = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) / det;
        inv[2][1] = -(jac[0][0] * jac[1][2] - jac[1][0] * jac[0][2]) / det;
        inv[0][2] = (jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1]) / det;
        inv[1][2] = -(jac[0][0] * jac[2][1] - jac[2][0] * jac[0][1]) / det;
        inv[2][2] = (jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1]) / det;
        break;
    }
    return std::make_pair(det, inv);
}

std::vector<std::vector<double> > solver::get_shape(type s, type e, type o, std::vector<double> pos) {
    std::vector<std::vector<double> > shp;
    switch (s) {
    case HGRAD:
        switch (e) {
        case EDGE:
            if (o == EYE) {
                shp = { { 1-pos[0] }, { pos[0]}};
                if(mdl.frm.p > 1)
                    shp.push_back( { 4*(1-pos[0])*pos[0]});
                if(mdl.frm.p > 2)
                    shp.push_back( { (1-2*pos[0])*(1-pos[0])*pos[0]});
            }
            break;
        case TRIA:
            switch(o) {
            case EYE:
                shp = { { 1-pos[0]-pos[1] }, { pos[0] }, { pos[1]}};
                if(mdl.frm.p>1)
                    shp.push_back( { { 4*pos[0]*pos[1] }, { 4*pos[1]*(1-pos[0]-pos[1]) }, { 4*pos[0]*(1-pos[0]-pos[1])}});
                if(mdl.frm.p>2)
                    shp.push_back( { { pos[0]*pos[1]*(pos[0]-pos[1]) },
                    { 	(1-pos[0]-pos[1])*pos[1]*(1-pos[0]-2*pos[1]) },
                    { 	(1-pos[0]-pos[1])*pos[0]*(1-2*pos[0]-pos[1]) },
                    { 	(1-pos[0]-pos[1])*pos[0]*pos[1]}
                });
                break;
            case GRAD:
                shp = { { -1,-1 }, { 1, 0 }, { 0, 1}};
                if(mdl.frm.p > 1) {
                    shp.push_back( { 4*pos[1], 4*pos[0]});
                    shp.push_back( { -4*pos[1], 4*(1-pos[0]-2*pos[1])});
                    shp.push_back( { 4*(1-2*pos[0]-pos[1]), -4*pos[0]});
                }
                break;
            }
            break;
        case TETRA:
            switch(o) {
            case EYE:
                shp = { { 1-pos[0]-pos[1]-pos[2] }, { pos[0] }, { pos[1] }, { pos[2]}};
                break;
            case GRAD:
                shp = { { -1,-1,-1 }, { 1,0,0 }, { 0,1,0 }, { 0,0,1}};
                break;
            }
            break;
        }
        break;
    case HCURL:
        switch(e) {
        case EDGE:
            break;
        case TRIA:
            break;
        case TETRA:
            break;
        }
        break;
    case HDIV:
        switch(e) {
        case EDGE:
            break;
        case TRIA:
            break;
        case TETRA:
            break;
        }
        break;
    }
    return shp;
}

std::vector<std::vector<double> > solver::get_stiff_mat(type s, type e, unsigned int i,
        std::pair<double, std::vector<std::vector<double> > >& jac) {
    std::vector<std::vector<double> > mat(i, std::vector<double>(i, 0));
    switch (s) {
    case HGRAD:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            if (mdl.frm.p == 1) {
                for (unsigned int iq = 0; iq < wq2_strang1.size(); iq++) {
                    auto shp = get_shape(HGRAD, TRIA, GRAD, xq2_strang1[iq]) * jac.second;
                    mat += (shp ^ shp) * wq2_strang1[iq];
                }
                mat *= jac.first;
            }
            break;
        case TETRA:
            if (mdl.frm.p == 1) {
                for (unsigned int iq = 0; iq < wq3_keast1.size(); iq++) {
                    auto shp = get_shape(HGRAD, TETRA, GRAD, xq3_keast1[iq]) * jac.second;
                    mat += (shp ^ shp) * wq3_keast1[iq];
                }
                mat *= jac.first;
            }
            break;
        }
        break;
    case HCURL:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            if (mdl.frm.p == 1) {
                for (unsigned int iq = 0; iq < wq2_strang5.size(); iq++) {
                    auto shp = get_shape(HCURL, TRIA, CURL, xq2_strang5[iq]) * jac.second;
                    mat += (shp ^ shp) * wq2_strang5[iq];
                }
                mat *= jac.first;
            }
            break;
        case TETRA:
            break;
        }
        break;
    case HDIV:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            break;
        case TETRA:
            break;
        }
        break;
    }
    return mat;
}
std::vector<std::vector<double> > solver::get_mass_mat(type s, type e, unsigned int i,
        std::pair<double, std::vector<std::vector<double> > >& jac) {
    std::vector<std::vector<double> > mat(i, std::vector<double>(i, 0));
    switch (s) {
    case HGRAD:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            break;
        case TETRA:
            break;
        }
        break;
    case HCURL:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            if (mdl.frm.p == 1) {
                for (unsigned int iq = 0; iq < wq2_strang5.size(); iq++) {
                    auto shp = get_shape(HCURL, TRIA, EYE, xq2_strang5[iq]) * jac.second;
                    mat += (shp ^ shp) * wq2_strang5[iq];
                }
                mat *= jac.first;
            }
            break;
            break;
        case TETRA:
            break;
        }
        break;
    case HDIV:
        switch (e) {
        case EDGE:
            break;
        case TRIA:
            break;
        case TETRA:
            break;
        }
        break;
    }
    return mat;
}

/// CUSTOM FORMULATION SOLVERS
void solver::analyze_em_e_fd() {
    sp_mat<double> matrix;
    sp_mat<double> rhs;
    double freq = mdl.frm.freq.range[0];
    double wavnbr = 2*phys_const::pi*freq/phys_const::c0;
    double wavnbr_2 = wavnbr*wavnbr;
    std::cout << "High frequency E-field formulation system assembly\n";
    std::vector<std::vector<bool> > dir(3);
    dir[HGRAD].assign(dof_num.hgrad, false);
    dir[HCURL].assign(dof_num.hcurl, false);
    for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
        if (strcmp(mdl.frm.bcs[i].type.data(), "PerfectE") == 0) {
            std::cout << mdl.frm.bcs[i].name << "\n";
            for (size_t j = 0; j < mdl.frm.bcs[i].faces.size(); j++) {
                auto dof = get_dofs(HCURL, TRIA, mdl.frm.bcs[i].faces[j]);
                for (unsigned int k = 0; k < dof.size(); k++) {
                    dir[HCURL][dof[k]] = true;
                }
                dof = get_dofs(HGRAD, TRIA, mdl.frm.bcs[i].faces[j]);
                for (unsigned int k = 0; k < dof.size(); k++) {
                    dir[HGRAD][dof[k]] = true;
                }
            }
        }
    }
    for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
        if (strcmp(mdl.frm.bcs[i].type.data(), "WavePort") == 0) {
            std::cout << mdl.frm.bcs[i].name << "\n";
            std::vector<size_t> h_curl_dofs, h_grad_dofs;
            for (size_t j = 0; j < mdl.frm.bcs[i].faces.size(); j++) {
                auto h_curl_dof = get_dofs(HCURL, TRIA, mdl.frm.bcs[i].faces[j]);
                auto h_grad_dof = get_dofs(HGRAD, TRIA, mdl.frm.bcs[i].faces[j]);
                h_curl_dofs.insert(std::end(h_curl_dofs), std::begin(h_curl_dof), std::end(h_curl_dof));
                h_grad_dofs.insert(std::end(h_grad_dofs), std::begin(h_grad_dof), std::end(h_grad_dof));
            }
            std::sort(h_curl_dofs.begin(), h_curl_dofs.end());
            std::sort(h_grad_dofs.begin(), h_grad_dofs.end());
            auto last = std::unique(h_curl_dofs.begin(), h_curl_dofs.end());
            h_curl_dofs.erase(last, h_curl_dofs.end());
            last = std::unique(h_grad_dofs.begin(), h_grad_dofs.end());
            h_grad_dofs.erase(last, h_grad_dofs.end());
            mdl.frm.bcs[i].mode_dof_map[HGRAD] = h_grad_dofs;
            mdl.frm.bcs[i].mode_dof_map[HCURL] = h_curl_dofs;
            double max_epsr = DBL_MIN;
            double max_mur = DBL_MIN;
            sp_mat<double> A, B;
            for (size_t j = 0; j < mdl.frm.bcs[i].faces.size(); j++) {
                auto id = mdl.frm.bcs[i].faces[j];
                double mur = mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[id]]].mur;
                double epsr = mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[id]]].epsr;
                max_epsr = std::max(max_epsr, epsr);
                max_mur = std::max(max_mur, mur);
                auto dof = get_dofs(HCURL, TRIA, id);
                auto jac = get_jac(mdl.msh.fac_geo(id));
                auto stiff = get_stiff_mat(HGRAD, TRIA, dof.size(), jac);
                auto mass = get_mass_mat(HGRAD, TRIA, dof.size(), jac);
                for (unsigned int i = 0; i < dof.size(); i++)
                    for (unsigned int j = 0; j < dof.size(); j++)
                        A(dof[i], dof[j]) += stiff[i][j]/mur - wavnbr_2*mass[i][j]*epsr ;
            }
            A.save("test.txt");
        }
    }
    /*
    if (mdl.msh.n_tetras == 0) {
        for (size_t id = 0; id < mdl.msh.n_faces; id++) {
            auto dof = get_dofs(HCURL, TRIA, id);
            auto jac = get_jac(mdl.msh.fac_geo(id));
            auto stiff = get_stiff_mat(HCURL, TRIA, dof.size(), jac)
                         / mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[id]]].mur;
            auto mass = get_mass_mat(HCURL, TRIA, dof.size(), jac)
                         * mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[id]]].epsr;
            for (unsigned int i = 0; i < dof.size(); i++)
                for (unsigned int j = 0; j < dof.size(); j++)
                    matrix(dof[i], dof[j]) += stiff[i][j] + ;
        }
        dir.assign(matrix.size(), false);
        for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
            if (strcmp(mdl.frm.bcs[i].type.data(), "Voltage") == 0) {
                for (size_t j = 0; j < mdl.frm.bcs[i].edges.size(); j++) {
                    auto dof = get_dofs(HGRAD, EDGE, mdl.frm.bcs[i].edges[j]);
                    for (unsigned int k = 0; k < dof.size(); k++) {
                        rhs(0, dof[k]) = mdl.frm.bcs[i].voltage;
                        dir[dof[k]] = true;
                    }
                }
            }
        }
    } else {
        for (size_t id = 0; id < mdl.msh.n_tetras; id++) {
            auto dof = get_dofs(HGRAD, TETRA, id);
            auto jac = get_jac(mdl.msh.tet_geo(id));
            auto stiff = get_stiff_mat(HGRAD, TETRA, dof.size(), jac)
                         * mdl.frm.mtrls[mtrl_map[mdl.msh.tet_lab[id]]].epsr;
            for (unsigned int i = 0; i < dof.size(); i++)
                for (unsigned int j = 0; j < dof.size(); j++)
                    matrix(dof[i], dof[j]) += stiff[i][j];
        }
        dir.assign(matrix.size(), false);
        for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
            if (strcmp(mdl.frm.bcs[i].type.data(), "Voltage") == 0) {
                for (size_t j = 0; j < mdl.frm.bcs[i].faces.size(); j++) {
                    auto dof = get_dofs(HGRAD, TRIA, mdl.frm.bcs[i].faces[j]);
                    for (unsigned int k = 0; k < dof.size(); k++) {
                        rhs(0, dof[k]) = mdl.frm.bcs[i].voltage;
                        dir[dof[k]] = true;
                    }
                }
            }
        }
    }
    if (rhs.size() == 0)
        throw std::string("Excitation not valid");
    std::cout << "Applying Dirichlet boundary conditions\n";
    sp_vec<double> sp_non_dir_val = matrix * rhs[0];
    for (size_t i = 0; i < dir.size(); i++)
        if (dir[i])
            matrix[i].clear();
    for (typename sp_mat<double>::iterator itr = matrix.begin();
            itr != matrix.end(); itr++)
        if (itr->size() > 0) {
            for (typename sp_vec<double>::iterator itc = itr->begin();
                    itc != itr->end(); itc++)
                if (dir[itc->idx])
                    itc->val = 0.0;
        }
    for (size_t i = 0; i < dir.size(); i++)
        if (dir[i])
            matrix(i, i) = 1.0;
    for (typename sp_vec<double>::iterator itc = sp_non_dir_val.begin();
            itc != sp_non_dir_val.end(); itc++)
        if (!dir[itc->idx])
            rhs(0, itc->idx) -= itc->val;
    std::cout << "Direct solver solution\n";
    solve_zmumps(matrix, rhs);
    */
}

void solver::analyze_e_v_stat() {
    sp_mat<double> matrix;
    sp_mat<double> rhs;
    std::vector<bool> dir;
    std::cout << "Electrostatic formulation system assembly\n";
    if (mdl.msh.n_tetras == 0) {
        for (size_t id = 0; id < mdl.msh.n_faces; id++) {
            auto dof = get_dofs(HGRAD, TRIA, id);
            auto jac = get_jac(mdl.msh.fac_geo(id));
            auto stiff = get_stiff_mat(HGRAD, TRIA, dof.size(), jac)
                         * mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[id]]].epsr;
            for (unsigned int i = 0; i < dof.size(); i++)
                for (unsigned int j = 0; j < dof.size(); j++)
                    matrix(dof[i], dof[j]) += stiff[i][j];
        }
        dir.assign(matrix.size(), false);
        for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
            if (strcmp(mdl.frm.bcs[i].type.data(), "Voltage") == 0) {
                for (size_t j = 0; j < mdl.frm.bcs[i].edges.size(); j++) {
                    auto dof = get_dofs(HGRAD, EDGE, mdl.frm.bcs[i].edges[j]);
                    for (unsigned int k = 0; k < dof.size(); k++) {
                        rhs(0, dof[k]) = mdl.frm.bcs[i].voltage;
                        dir[dof[k]] = true;
                    }
                }
            }
        }
    } else {
        for (size_t id = 0; id < mdl.msh.n_tetras; id++) {
            auto dof = get_dofs(HGRAD, TETRA, id);
            auto jac = get_jac(mdl.msh.tet_geo(id));
            auto stiff = get_stiff_mat(HGRAD, TETRA, dof.size(), jac)
                         * mdl.frm.mtrls[mtrl_map[mdl.msh.tet_lab[id]]].epsr;
            for (unsigned int i = 0; i < dof.size(); i++)
                for (unsigned int j = 0; j < dof.size(); j++)
                    matrix(dof[i], dof[j]) += stiff[i][j];
        }
        dir.assign(matrix.size(), false);
        for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
            if (strcmp(mdl.frm.bcs[i].type.data(), "Voltage") == 0) {
                for (size_t j = 0; j < mdl.frm.bcs[i].faces.size(); j++) {
                    auto dof = get_dofs(HGRAD, TRIA, mdl.frm.bcs[i].faces[j]);
                    for (unsigned int k = 0; k < dof.size(); k++) {
                        rhs(0, dof[k]) = mdl.frm.bcs[i].voltage;
                        dir[dof[k]] = true;
                    }
                }
            }
        }
    }
    if (rhs.size() == 0)
        throw std::string("Excitation not valid");
    std::cout << "Applying Dirichlet boundary conditions\n";
    sp_vec<double> sp_non_dir_val = matrix * rhs[0];
    for (size_t i = 0; i < dir.size(); i++)
        if (dir[i])
            matrix[i].clear();
    for (typename sp_mat<double>::iterator itr = matrix.begin();
            itr != matrix.end(); itr++)
        if (itr->size() > 0) {
            for (typename sp_vec<double>::iterator itc = itr->begin();
                    itc != itr->end(); itc++)
                if (dir[itc->idx])
                    itc->val = 0.0;
        }
    for (size_t i = 0; i < dir.size(); i++)
        if (dir[i])
            matrix(i, i) = 1.0;
    for (typename sp_vec<double>::iterator itc = sp_non_dir_val.begin();
            itc != sp_non_dir_val.end(); itc++)
        if (!dir[itc->idx])
            rhs(0, itc->idx) -= itc->val;
    std::cout << "Direct solver solution\n";
    solve_dmumps(matrix, rhs);
}

void solver::analyze_h_a_stat() {
    sp_mat<double> matrix;
    sp_mat<double> rhs;
    std::vector<bool> dir;
    std::cout << "Magnetostatic formulation system assembly\n";
    if (mdl.msh.n_tetras == 0) {
        for (size_t id = 0; id < mdl.msh.n_faces; id++) {
            auto dof = get_dofs(HGRAD, TRIA, id);
            auto jac = get_jac(mdl.msh.fac_geo(id));
            auto stiff = get_stiff_mat(HGRAD, TRIA, dof.size(), jac)
                         * mdl.frm.mtrls[mtrl_map[mdl.msh.fac_lab[id]]].epsr;
            for (unsigned int i = 0; i < dof.size(); i++)
                for (unsigned int j = 0; j < dof.size(); j++)
                    matrix(dof[i], dof[j]) += stiff[i][j];
        }
        dir.assign(matrix.size(), false);
        for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
            if (strcmp(mdl.frm.bcs[i].type.data(), "Current") == 0) {
                for (size_t j = 0; j < mdl.frm.bcs[i].edges.size(); j++) {
                    auto dof = get_dofs(HGRAD, EDGE, mdl.frm.bcs[i].edges[j]);
                    for (unsigned int k = 0; k < dof.size(); k++) {
                        rhs(0, dof[k]) = mdl.frm.bcs[i].current;
                        dir[dof[k]] = true;
                    }
                }
            }
        }
    } else {
        for (size_t id = 0; id < mdl.msh.n_tetras; id++) {
            auto dof = get_dofs(HGRAD, TETRA, id);
            auto jac = get_jac(mdl.msh.tet_geo(id));
            auto stiff = get_stiff_mat(HGRAD, TETRA, dof.size(), jac)
                         * mdl.frm.mtrls[mtrl_map[mdl.msh.tet_lab[id]]].epsr;
            for (unsigned int i = 0; i < dof.size(); i++)
                for (unsigned int j = 0; j < dof.size(); j++)
                    matrix(dof[i], dof[j]) += stiff[i][j];
        }
        dir.assign(matrix.size(), false);
        for (size_t i = 0; i < mdl.frm.bcs.size(); i++) {
            if (strcmp(mdl.frm.bcs[i].type.data(), "Current") == 0) {
                for (size_t j = 0; j < mdl.frm.bcs[i].faces.size(); j++) {
                    auto dof = get_dofs(HGRAD, TRIA, mdl.frm.bcs[i].faces[j]);
                    for (unsigned int k = 0; k < dof.size(); k++) {
                        rhs(0, dof[k]) = mdl.frm.bcs[i].voltage;
                        dir[dof[k]] = true;
                    }
                }
            }
        }
    }
    if (rhs.size() == 0)
        throw std::string("Excitation not valid");
    std::cout << "Applying Dirichlet boundary conditions\n";
    sp_vec<double> sp_non_dir_val = matrix * rhs[0];
    for (size_t i = 0; i < dir.size(); i++)
        if (dir[i])
            matrix[i].clear();
    for (typename sp_mat<double>::iterator itr = matrix.begin();
            itr != matrix.end(); itr++)
        if (itr->size() > 0) {
            for (typename sp_vec<double>::iterator itc = itr->begin();
                    itc != itr->end(); itc++)
                if (dir[itc->idx])
                    itc->val = 0.0;
        }
    for (size_t i = 0; i < dir.size(); i++)
        if (dir[i])
            matrix(i, i) = 1.0;
    for (typename sp_vec<double>::iterator itc = sp_non_dir_val.begin();
            itc != sp_non_dir_val.end(); itc++)
        if (!dir[itc->idx])
            rhs(0, itc->idx) -= itc->val;
    std::cout << "Direct solver solution\n";
    solve_dmumps(matrix, rhs);
}

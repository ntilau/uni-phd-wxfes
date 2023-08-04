#ifndef SOLVER_H
#define SOLVER_H

#include "project.h"

#include <vector>
#include <map>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>

template<typename T>
class sp_node {
public:
    sp_node(size_t i) :
        idx(i), val(0) {
    }
    ~sp_node() {
    }
    bool operator==(sp_node<T> tmp) {
        return idx == tmp.idx;
    }
    bool operator==(const sp_node<T>& tmp) const {
        return idx == tmp.idx;
    }
    T val;
    size_t idx;
};

template<typename T>
class sp_vec: public std::vector<sp_node<T> > {
public:
    typedef typename std::vector<sp_node<T> >::iterator iterator;
    typedef typename std::vector<sp_node<T> >::const_iterator const_iterator;
    T& operator[](size_t i) {
        iterator it = std::find(this->begin(), this->end(), sp_node<T>(i));
        if (it == this->end()) {
            this->push_back(sp_node<T>(i));
            return this->back().val;
        }
        return it->val;
    }
    T operator[](size_t i) const {
        const_iterator it = std::find(this->begin(), this->end(),
                                      sp_node<T>(i));
        if (it == this->end())
            return T(0);
        else
            return it->val;
    }
    std::vector<T> full() {
        std::vector<T> tmp(this->dim(), 0);
        for (iterator it = this->begin(); it != this->end(); it++)
            tmp[it->idx] = it->val;
        return tmp;
    }
    size_t nnz() {
        size_t tmp = 0;
        for (iterator it = this->begin(); it != this->end(); it++)
            ++tmp;
        return tmp;
    }
    size_t dim() {
        size_t tmp = 0;
        for (iterator it = this->begin(); it != this->end(); it++)
            tmp = std::max(tmp, it->idx);
        return tmp + 1;
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& out, sp_vec<T> v) {
    out << "sparse vector\n";
    for (typename sp_vec<T>::iterator it = v.begin(); it != v.end(); it++) {
        out << "[" << it->idx << "]=" << it->val << "\n";
    }
    return out;
}

template<typename T>
class sp_mat: public std::vector<sp_vec<T> > {
public:
    typedef typename std::vector<sp_vec<T> >::iterator iterator;
    T& operator()(size_t i, size_t j) {
        if (i >= this->size()) {
            for (size_t row = this->size(); row <= i; row++)
                this->push_back(sp_vec<T>());
        }
        return this->operator[](i)[j];
    }
    sp_vec<T> operator*(const sp_vec<T>& vec) {
        sp_vec<T> tmp;
        for (typename sp_mat<T>::iterator itr = this->begin();
                itr != this->end(); itr++) {
            for (typename sp_vec<T>::iterator itc = itr->begin();
                    itc != itr->end(); itc++) {
                T val = itc->val * vec[itc->idx];
                if (val != T(0))
                    tmp[itr - this->begin()] += val;
            }
        }
        return tmp;
    }
    size_t nnz() {
        size_t tmp = 0;
        for (typename sp_mat<T>::iterator itr = this->begin();
                itr != this->end(); itr++)
            tmp += itr->nnz();
        return tmp;
    }
    void save(std::string name) {
        std::ofstream mat_file(name.data());
        for (typename sp_mat<T>::iterator itr = this->begin();
                itr != this->end(); itr++) {
            for (typename sp_vec<T>::iterator itc = itr->begin();
                    itc != itr->end(); itc++) {
                mat_file << itr - this->begin() + 1 << " " << itc->idx + 1
                         << " " << std::setprecision(16) << itc->val << "\n";
            }
        }
        mat_file.close();
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& out, sp_mat<T> v) {
    out << "sparse matrix\n";
    for (typename sp_mat<T>::iterator itr = v.begin(); itr != v.end(); itr++) {
        for (typename sp_vec<T>::iterator itc = itr->begin(); itc != itr->end();
                itc++) {
            out << "[" << itr - v.begin() << "][" << itc->idx << "]="
                << itc->val << "\n";
        }
    }
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream &out, const std::vector<T> &vec) {
    // out << "matrix " << mat.val.size() << "\n";
    out << "{ ";
    for (size_t i = 0; i < vec.size(); i++)
        out << "[" << i << "] " << vec[i] << " ";
    out << "}\n";
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream &out,
                         const std::vector<std::vector<T> > &mat) {
    // out << "matrix " << mat.val.size() << "\n";
    out << "{ ";
    for (size_t i = 0; i < mat.size(); i++)
        out << "[" << i << "] " << mat[i] << " ";
    out << "}\n";
    return out;
}

template<typename T>
std::vector<std::vector<T> > operator*(const std::vector<std::vector<T> >& l,
                                       const std::vector<std::vector<T> >& r) {
    unsigned int m = l.size(), n = r[0].size(), p = r.size();
    std::vector<std::vector<T> > tmp(m, std::vector<T>(n, 0));
    for (unsigned int i = 0; i < m; i++)
        for (unsigned int j = 0; j < n; j++)
            for (unsigned int k = 0; k < p; k++)
                tmp[i][j] += l[i][k] * r[k][j];
    return tmp;
}

template<typename T>
std::vector<std::vector<T> > operator^(const std::vector<std::vector<T> >& l,
                                       const std::vector<std::vector<T> >& r) {
    unsigned int m = l.size(), n = r.size(), p = r[0].size();
    std::vector<std::vector<T> > tmp(m, std::vector<T>(n, 0));
    for (unsigned int i = 0; i < m; i++)
        for (unsigned int j = 0; j < n; j++)
            for (unsigned int k = 0; k < p; k++)
                tmp[i][j] += l[i][k] * r[j][k];
    return tmp;
}

template<typename T>
std::vector<std::vector<T> >& operator+=(std::vector<std::vector<T> >& l,
        const std::vector<std::vector<T> >& r) {
    unsigned int m = l.size(), n = l[0].size();
    for (unsigned int i = 0; i < m; i++)
        for (unsigned int j = 0; j < n; j++)
            l[i][j] += r[i][j];
    return l;
}

template<typename T>
std::vector<std::vector<T> > operator*(const std::vector<std::vector<T> >& l,
                                       const T& r) {
    unsigned int m = l.size(), n = l[0].size();
    std::vector<std::vector<T> > tmp(l);
    for (unsigned int i = 0; i < m; i++)
        for (unsigned int j = 0; j < n; j++)
            tmp[i][j] *= r;
    return tmp;
}

template<typename T>
std::vector<std::vector<T> >& operator*=(std::vector<std::vector<T> >& l,
        const T& r) {
    unsigned int m = l.size(), n = l[0].size();
    for (unsigned int i = 0; i < m; i++)
        for (unsigned int j = 0; j < n; j++)
            l[i][j] *= r;
    return l;
}

extern "C" void
dnaupd_(int* ido, char* problem, int* n, char* whichArr, int* nev, double* tol,
        double* residArr, int* ncv, double* vMatrix, int* ldv, int* iparamArr,
        int* ipntrArr, double* workdArr, double* worklArr, int* lworkl,
        int* info);

extern "C" void
dneupd_(int* rvec, char* stringArr, int* selectArr, double* dArrReal,
        double* dArrImag, double* vMatrix, int* ldv, double* sigmaReal,
        double* sigmaImag, double* workev, char* bmat, int* n, char* whichArr,
        int* nev, double* tol, double* residArr, int* ncv, double* vMatrix1,
        int* ldv1, int* iparamArr, int* ipntrArr, double* workdArr,
        double* worklArr, int* lworkl, int* ierr);

extern "C" void
znaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
        std::complex<double>* resid, int* ncv, std::complex<double>* v,
        int* ldv, int* iparam, int* ipntr, std::complex<double>* workd,
        std::complex<double>* workl, int* lworkl, double* rwork, int* info);

extern "C" void
zneupd_(int* rvec, char* All, int* select, std::complex<double>* d,
        std::complex<double>* v, int* ldv, double* sigma,
        std::complex<double>* workev, char* bmat, int* n, char* which, int* nev,
        double* tol, std::complex<double>* resid, int* ncv,
        std::complex<double>* v1, int* ldv1, int* iparam, int* ipntr,
        std::complex<double>* workd, std::complex<double>* workl, int* lworkl,
        double* rwork, int* ierr);

class solver {
public:
    solver(mdl_core& mdl);
    virtual ~solver();
protected:
    enum type {
        HGRAD, HCURL, HDIV, EYE, GRAD, CURL, DIV, EDGE, TRIA, TETRA
    };
    void compute_label_map();
    void compute_dof_num();
    void solve_eig_zmumps(sp_mat<std::complex<double> >& a, sp_mat<std::complex<double> >& b,
                          size_t& n_t, size_t& n_z, double& kk, int& n_modes);
    void solve_dmumps(sp_mat<double>& matrix, sp_mat<double>& rhs);
    void solve_zmumps(sp_mat<std::complex<double> >& matrix,
                      sp_mat<std::complex<double> >& rhs);
    std::vector<size_t> get_dofs(type s, type e, size_t id);
    std::pair<double, std::vector<std::vector<double> > > get_jac(std::vector<std::vector<double> > geo);
    std::vector<std::vector<double> > get_shape(type s, type e, type o, std::vector<double> pos);
    std::vector<std::vector<double> > get_stiff_mat(type s, type e, unsigned int i, std::pair<double, std::vector<std::vector<double> > >& jac);
    std::vector<std::vector<double> > get_mass_mat(type s, type e, unsigned int i, std::pair<double, std::vector<std::vector<double> > >& jac);
    void analyze_em_e_eig();
    void analyze_em_e_fd();
    void analyze_e_v_stat();
    void analyze_h_a_stat();
private:
    mdl_core& mdl;
    struct DOFS {
        size_t hgrad = 0, hcurl = 0, hdiv = 0;
    } dof_num;
    std::map<int, unsigned int> mtrl_map, bc_map;
};

#endif // SOLVER_H

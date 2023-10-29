#ifndef MODEL_H
#define MODEL_H

#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <string>
#include <vector>
#include <list>
#include <complex>
#include <cstdlib>
#include <cfloat>
#include <climits>

class phys_const {
public:
    static constexpr double c0 = 299792458;
    static constexpr double eps0 = 8.854187817e-12;
    static constexpr double mu0 = 1.2566370614359e-6;
    static constexpr double pi = 3.1415926535898;
    static constexpr double z0 = 3.767303134749663e+02;
    static constexpr double freq = 9.4e+009;
};

class mdl_sld {
public:
    mdl_sld();
    ~mdl_sld();
    void clear();
    bool read_poly_file(std::string& name);
    bool read_stl_file(std::string& name);
    void write_prj_file(std::string& name);
    void read_prj_file(std::string& name);
    unsigned int check_dim();
    double get_geom_dim();
    double max_dimension = 0.0;
    std::vector<std::vector<double> > bounding_box;
    void get_bounding_info();
    /// poly based modeling primitives
    std::vector<std::vector<double> > nodes;
    std::vector<std::vector<size_t> > edges;
    std::vector<std::vector<int> > edges_marker;
    struct face_type {
        std::vector<std::vector<size_t> > polygons;
        std::vector<std::vector<double> > holes;
    };
    std::vector<face_type> faces;
    std::vector<std::vector<int> > faces_marker;
    std::vector<std::vector<double> > faces_normals;
    std::vector<std::vector<double> > holes;
    std::vector<std::vector<double> > regions;
    std::vector<int> regions_marker;
    int max_edges_marker = -SIZE_MAX,
        max_faces_marker = -SIZE_MAX,
        max_regions_marker = -SIZE_MAX;
    unsigned int dim = 0;
    std::string tetgen_command = "pqaAfee";
    std::string triangle_command = "pqaAe";
    std::vector<int> bc_markers;
    std::vector<int> mtrl_markers;
    std::string sld_name;
};

class mdl_bc {
public:
    mdl_bc();
    ~mdl_bc();
    int label;
    std::string name = "None";
    std::string type = "None";
    double power = 1.0;
    unsigned int num_modes = 0;
    bool tfe = true;
    std::vector<std::vector<size_t> > mode_dof_map = std::vector<std::vector<size_t> >(3); // HGRAD, HCURL, HDIV
    std::vector<std::complex<double> > mode_beta;
    std::vector<std::vector<std::complex<double> > > mode_eig_vec;
    std::vector<std::vector<std::complex<double> > > mode_eig_vec_f;
    std::complex<double> lumped_impedance = std::complex<double>(50.0, 0.0);
    double R = 50.0, L = 0.0, C = 0.0;
    std::complex<double> surf_impedance = std::complex<double>(0.0, 0.0);
    bool inc = false;
    std::vector<double> inc_E = { 1.0, 0.0, 0.0 };
    std::vector<double> inc_k = { 0.0, 0.0, 1.0 };
    double n_theta = 101, n_phi = 201;
    double voltage = 1.0;
    double current = 1.0;
    std::vector<size_t> faces;
    std::vector<size_t> edges;
};

class mdl_mtrl {
public:
    mdl_mtrl();
    ~mdl_mtrl();
    void upd_mtrl(double& freq);
    double calc_epsr2(double& freq);
    void upd_mtrl();
    double epsr = 1.0;
    double epsr2 = 0.0;
    double mur = 1.0;
    double kr = 0.0;
    double sigma = 0.0;
    double etaSigma = 0.0;
    double tand = 0.0;
    double kerr = 0.0;
    std::string name = "None";
    std::string type = "Vacuum";
    int label = 0;
    std::vector<size_t> tetras;
    std::vector<size_t> faces;
};

class mdl_msh {
public:
    mdl_msh();
    mdl_msh(mdl_msh*); // copy constructor
    ~mdl_msh();
    std::vector<std::string> mesh_type = { "TRIA", "TETRA" };
    std::string type = "TETRA";
    // methods
    double get_geom_dim();
    double max_dimension = 0.0;
    std::vector<std::vector<double> > bounding_box;
    void get_bounding_info();
    void write_prj_file(std::string& name);
    void read_prj_file(std::string& name);
    void read_tetgen_files(std::string& name);
    void read_triangle_files(std::string& name);
    void save_vtk_mesh(std::string);
    void regularize_mesh();
    void refine_homogeneous();
    void get_mesh_statistics();
    void clear();
    //
    std::vector<std::vector<double> > tet_geo(size_t);
    std::vector<std::vector<double> > fac_geo(size_t);
    std::vector<std::vector<double> > fac_geo2(size_t);
    std::vector<std::vector<double> > edg_geo(size_t);
    std::vector<double> int_node(size_t);
    std::vector<double> int_node(size_t, size_t&);
    // members
    std::vector<std::vector<size_t> > tet_nodes;
    std::vector<std::vector<size_t> > tet_edges;
    std::vector<std::vector<size_t> > tet_faces;
    std::vector<std::vector<size_t> > fac_nodes;
    std::vector<std::vector<size_t> > fac_edges;
    std::vector<std::vector<size_t> > edg_nodes;
    std::vector<std::vector<double> > nod_pos;
    std::vector<int> edg_lab;
    std::vector<int> fac_lab;
    std::vector<int> tet_lab;
    std::vector<std::vector<size_t> > fac_adj_tet;
    std::vector<std::vector<size_t> > edg_adj_fac;
    std::vector<std::vector<size_t> > dom_tetras;
    std::vector<std::vector<size_t> > dom_faces;
    size_t n_nodes, n_edges, n_faces, n_tetras, n_domains;
    int max_edg_marker = -SIZE_MAX,
        max_fac_marker = -SIZE_MAX,
        max_tet_marker = -SIZE_MAX;
};

class mdl_src {
public:
    mdl_src();
    mdl_src(mdl_src*); // copy constructor
    ~mdl_src();
    struct dipole {
        std::vector<std::complex<double> > amplitude;
        std::vector<double> position;
        std::vector<double> direction;
        double length;
    };
    std::vector<dipole> currents;
    std::string name = "None";
    std::string type = "None";
};

class mdl_frm {
public:
    std::vector<std::list<std::string> > bc_type = {
        { "None", "PerfectE", "PerfectH", "Radiation", "WavePort", "Impedance", "LumpedPort", "LumpedRLC"},
        { "None", "Voltage" },
        { "None", "Current", "Voltage", "Insulation", "Skin" }
    };
    std::vector<std::list<std::string> > frm_type = {
        { "EM_E_FD"  },
        { "E_V_STAT" },
        { "H_A_STAT" },
        { "EM_PO"    }
    };
    std::string type = "EM_E_FD";
    mdl_frm();
    ~mdl_frm();
    void clear();
    void update_msh_info(mdl_msh& msh);
    std::vector<mdl_bc> bcs;
    std::vector<mdl_mtrl> mtrls;
    std::vector<mdl_src> srcs;
    void write_prj_file(std::string& name);
    void read_prj_file(std::string& name);
    struct freq_type {
        std::vector<double> range = { 1e9, 1e9 };
        unsigned int nbr = 1;
    } freq;
    unsigned int niter = 100; // max number of iterations of iterative solver
    double toll = 1e-6; // tollerance for iterative solver
    double relax = 1.0;
    unsigned int h = 0; // homogeneous refinement
    unsigned int p = 1; // polynomial order
    /// post processing
    std::vector<std::vector<double> > sol_real;
    std::vector<std::vector<std::complex<double> > > sol_cmplx;
};

class mdl_core {
public:
    mdl_frm frm;
    mdl_sld sld;
    mdl_msh msh;
    void wrap_hfss(std::string& data_path, std::string& name);
    void wrap_aedt(std::string& data_path, std::string& name);
    void create_tri_mesh();
    void clear() {
        sld.clear();
        frm.clear();
        msh.clear();
    }
protected:
    void removeCharsFromString(std::string &str, const char* charsToRemove) {
        for (unsigned int i = 0; i < strlen(charsToRemove); ++i) {
            str.erase(std::remove(str.begin(), str.end(), charsToRemove[i]),
                      str.end());
        }
    }
    char find_SI_factor(std::string &str) {
        for (unsigned int i = 0; i < strlen(SI_chars); ++i) {
            size_t found = str.find(SI_chars[i]);
            if (found <= str.size())
                return str[found];
        }
        return 0;
    }
    double set_factor(char fact) {
        if (fact == 'm')
            return 1e-3;
        else if (fact == 'u')
            return 1e-6;
        else if (fact == 'n')
            return 1e-9;
        else if (fact == 'p')
            return 1e-12;
        else
            return 1;
    }
    const char* SI_chars = "munp";
};

#endif // MODEL_H

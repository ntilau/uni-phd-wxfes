#include "model.h"
#include <map>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <cstring>

mdl_sld::mdl_sld() {
}

mdl_sld::~mdl_sld() {
}

double mdl_sld::get_geom_dim() {
    double dim = 0;
    for (unsigned int i = 0; i < nodes.size(); i++) {
        dim = std::max(dim, std::abs(nodes[i][0]));
        dim = std::max(dim, std::abs(nodes[i][1]));
        dim = std::max(dim, std::abs(nodes[i][2]));
    }
    dim *= 1e-1;
    for (int i = -15; i < 15; i++) {
        if (dim < pow(10, double(i)))
            return pow(10, double(i));
    }
    return 0;
}

bool mdl_sld::read_stl_file(std::string& name) {
    clear();
    std::map<std::tuple<double, double, double>, size_t> nod_map;
    bool is_binary = false;
    std::string tmp_str, tmp_str_sld;
    std::string line;
    std::vector<double> tmp_pnt(3);
    std::vector<size_t> tmp_tri(3);
    size_t nod_cnt = 0, fac_cnt = 0;
    std::ifstream stl_file(std::string(name + ".stl").data());
    if (stl_file.is_open()) {
        getline(stl_file, line);
        std::istringstream iss;
        iss.str(line);
        iss >> tmp_str_sld;
        if (tmp_str_sld == "solid") {
            iss >> sld_name;
            std::cout << "Solid name: " << sld_name << "\n";
            dim = 3;
            do {
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // facet
                if(tmp_str == "endsolid") {
                    getline(stl_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> tmp_str; // facet
                    if(tmp_str == "endsolid")
                        break;
                    else if(stl_file.eof())
                        break;
                }
                iss >> tmp_str; // normal
                iss >> tmp_pnt[0]; // ni
                iss >> tmp_pnt[1]; // nj
                iss >> tmp_pnt[2]; // nk
                faces_normals.push_back(tmp_pnt);
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // outer
                iss >> tmp_str; // loop
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // vertex
                iss >> tmp_pnt[0]; // v1x
                iss >> tmp_pnt[1]; // v1y
                iss >> tmp_pnt[2]; // v1z
                if (nod_map.find(std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])) == nod_map.end()) {
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])] =  nod_cnt++;
                }
                tmp_tri[0] = nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])];
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // vertex
                iss >> tmp_pnt[0]; // v2x
                iss >> tmp_pnt[1]; // v2y
                iss >> tmp_pnt[2]; // v2z
                if (nod_map.find(std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])) == nod_map.end()) {
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])] = nod_cnt++;
                }
                tmp_tri[1] = nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])];
                getline(stl_file, line);
                iss.clear();
                iss.str(line);
                iss >> tmp_str; // vertex
                iss >> tmp_pnt[0]; // v3x
                iss >> tmp_pnt[1]; // v3y
                iss >> tmp_pnt[2]; // v3z
                if (nod_map.find(std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])) == nod_map.end()) {
                    nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])] = nod_cnt++;
                }
                tmp_tri[2] = nod_map[std::make_tuple(tmp_pnt[0], tmp_pnt[1], tmp_pnt[2])];
                face_type fac_ply;
                fac_ply.polygons.push_back(tmp_tri);
                faces.push_back(fac_ply);
                faces_marker.push_back(std::vector<int>(1,1));
                fac_cnt++;
                getline(stl_file, line); // endloop
                getline(stl_file, line); // endfacet
            } while (!stl_file.eof());
            stl_file.close();
            std::cout << "Nodes = " << nod_cnt << "\n";
            std::cout << "Faces = " << fac_cnt << "\n";
            nodes.resize(nod_cnt);
            for (std::map<std::tuple<double, double, double>, size_t>::iterator it =
                        nod_map.begin(); it != nod_map.end(); it++) {
                std::vector<double> node(3);
                node[0] = std::get < 0 > (it->first);
                node[1] = std::get < 1 > (it->first);
                node[2] = std::get < 2 > (it->first);
                nodes[it->second] = node;
            }
            max_faces_marker = 1;
            bc_markers.push_back(1);
        } else {
            is_binary = true;
        }
    } else {
        std::cout << name + ".stl not found\n";
        return false;
    }
    return true;
}

bool mdl_sld::read_poly_file(std::string& name) {
    clear();
    size_t tmp_int, n_pts;
    double tmp_dbl;
    std::string tmp_str;
    std::string line;
    std::ifstream poly_file(std::string(name + ".poly").data());
    if (poly_file.is_open()) {
        {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            nodes.resize(tmp_int);
            iss >> n_pts;
            for (size_t i = 0; i < nodes.size(); i++) {
                do {
                    getline(poly_file, line);
                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                std::istringstream iss(line);
                iss >> tmp_int;
                iss >> tmp_dbl;
                nodes[i].push_back(tmp_dbl);
                iss >> tmp_dbl;
                nodes[i].push_back(tmp_dbl);
                iss >> tmp_dbl;
                if (n_pts < 3)
                    nodes[i].push_back(0.0);
                else
                    nodes[i].push_back(tmp_dbl);
            }
        }
        dim = check_dim();
        std::cout << "Model is " << dim << "-dimensional\n";
        if (dim == 3) {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            faces.resize(tmp_int);
            faces_marker.resize(tmp_int);
            unsigned int n_faces_marker;
            iss >> n_faces_marker;
            for (size_t i = 0; i < faces.size(); i++) {
                int polygons = 0, hole = 0, bmark = 0;
                do {
                    getline(poly_file, line);
                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                std::istringstream iss(line);
                iss >> polygons;
                iss >> hole;
                iss >> bmark;
                max_faces_marker = std::max(max_faces_marker, bmark);
                if (std::find(bc_markers.begin(), bc_markers.end(), bmark)
                        == bc_markers.end())
                    bc_markers.push_back(bmark);
                faces_marker[i].push_back(bmark);
                faces[i].polygons.resize(polygons);
                if (hole > 0)
                    faces[i].holes.resize(hole);
                for (size_t j = 0; j < polygons; j++) {
                    unsigned int nodes_nbr = 0;
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> nodes_nbr;
                    for (unsigned int k = 0; k < nodes_nbr; k++) {
                        iss >> tmp_int;
                        if (faces[i].polygons[j].size() > 0) {
                            if (faces[i].polygons[j].back() == (tmp_int - 1)) {
                                do {
                                    getline(poly_file, line);
                                } while (line.compare(0, 1, "#") == 0
                                         || line.size() == 0);
                                iss.clear();
                                iss.str("");
                                iss.str(line);
                                iss >> tmp_int;
                            }
                        }
                        faces[i].polygons[j].push_back(tmp_int - 1);
                    }
                }
                for (size_t j = 0; j < hole; j++) {
                    unsigned int nodes_nbr = 0;
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> tmp_int;
                    iss >> tmp_dbl;
                    faces[i].holes[j].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    faces[i].holes[j].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    faces[i].holes[j].push_back(tmp_dbl);
                }
            }
        } else if (dim == 2) {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            edges.resize(tmp_int);
            edges_marker.resize(tmp_int);
            iss >> tmp_int;
            size_t n_edges_marker = tmp_int;
            for (size_t i = 0; i < edges.size(); i++) {
                do {
                    getline(poly_file, line);
                } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                std::istringstream iss(line);
                iss >> tmp_int;
                iss >> tmp_int;
                edges[i].push_back(tmp_int - 1);
                iss >> tmp_int;
                edges[i].push_back(tmp_int - 1);
                for (unsigned int j = 0; j < n_edges_marker; j++) {
                    iss >> tmp_int;
                    int bmark = tmp_int;
                    max_edges_marker = std::max(max_edges_marker, bmark);
                    edges_marker[i].push_back(bmark);
                    if (std::find(bc_markers.begin(), bc_markers.end(), bmark)
                            == bc_markers.end())
                        bc_markers.push_back(bmark);
                }
            }

        }
        // importing holes
        {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            if (tmp_int > 0) {
                holes.resize(tmp_int);
                for (size_t i = 0; i < holes.size(); i++) {
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> tmp_int;
                    iss >> tmp_dbl;
                    holes[i].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    holes[i].push_back(tmp_dbl);
                    if (dim == 2)
                        holes[i].push_back(0.0);
                }
            }
        }
        // importing regions
        {
            do {
                getline(poly_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            std::istringstream iss(line);
            iss >> tmp_int;
            if (tmp_int > 0) {
                regions.resize(tmp_int);
                regions_marker.resize(tmp_int);
                for (size_t i = 0; i < regions.size(); i++) {
                    do {
                        getline(poly_file, line);
                    } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
                    std::istringstream iss(line);
                    iss >> tmp_int;
                    iss >> tmp_dbl;
                    regions[i].push_back(tmp_dbl);
                    iss >> tmp_dbl;
                    regions[i].push_back(tmp_dbl);
                    if (dim == 3) {
                        iss >> tmp_dbl;
                        regions[i].push_back(tmp_dbl);
                    } else if (dim == 2)
                        regions[i].push_back(0.0);
                    iss >> regions_marker[i];
                    max_regions_marker = std::max(max_regions_marker,
                                                  regions_marker[i]);
                    if (std::find(mtrl_markers.begin(), mtrl_markers.end(),
                                  regions_marker[i]) == mtrl_markers.end())
                        mtrl_markers.push_back(regions_marker[i]);
                }
            }
        }
        std::sort(bc_markers.begin(), bc_markers.end());
        std::sort(mtrl_markers.begin(), mtrl_markers.end());
        return true;
    } else {
        std::cout << name + ".poly not found\n";
        return false;
    }
}

void mdl_sld::clear() {
    std::vector<std::vector<double> >().swap(nodes);
    std::vector<std::vector<size_t> >().swap(edges);
    std::vector<std::vector<int> >().swap(edges_marker);
    std::vector<face_type>().swap(faces);
    std::vector<std::vector<int> >().swap(faces_marker);
    std::vector<std::vector<double> >().swap(holes);
    std::vector<std::vector<double> >().swap(regions);
    std::vector<int>().swap(regions_marker);
    std::vector<std::vector<double> >().swap(bounding_box);
    std::vector<int>().swap(bc_markers);
    std::vector<int>().swap(mtrl_markers);
    max_dimension = 0.0;
    dim = 0;
    max_edges_marker = -INT_MAX;
    max_faces_marker = -INT_MAX;
    max_regions_marker = -INT_MAX;
}

void mdl_sld::get_bounding_info() {
    bounding_box.resize(2);
    bounding_box[0].assign(3, DBL_MAX);
    bounding_box[1].assign(3, -DBL_MAX);
    for (size_t i = 0; i < nodes.size(); i++) {
        bounding_box[0][0] = std::min(bounding_box[0][0], nodes[i][0]);
        bounding_box[0][1] = std::min(bounding_box[0][1], nodes[i][1]);
        bounding_box[0][2] = std::min(bounding_box[0][2], nodes[i][2]);
        bounding_box[1][0] = std::max(bounding_box[1][0], nodes[i][0]);
        bounding_box[1][1] = std::max(bounding_box[1][1], nodes[i][1]);
        bounding_box[1][2] = std::max(bounding_box[1][2], nodes[i][2]);
    }
    double x_dim, y_dim, z_dim;
    x_dim = std::abs(bounding_box[1][0] - bounding_box[0][0]);
    y_dim = std::abs(bounding_box[1][1] - bounding_box[0][1]);
    z_dim = std::abs(bounding_box[1][2] - bounding_box[0][2]);
    max_dimension = std::max(max_dimension, x_dim);
    max_dimension = std::max(max_dimension, y_dim);
    max_dimension = std::max(max_dimension, z_dim);
}

unsigned int mdl_sld::check_dim() {
    double x, y, z;
    if (nodes.size() >= 1) {
        x = nodes[0][0];
        y = nodes[0][1];
        z = nodes[0][2];
    }
    double x_norm = 0, y_norm = 0, z_norm = 0;
    for (size_t i = 0; i < nodes.size(); i++) {
        x_norm += std::abs(nodes[i][0] * nodes[i][0]);
        y_norm += std::abs(nodes[i][1] * nodes[i][1]);
        z_norm += std::abs(nodes[i][2] * nodes[i][2]);
    }
    x_norm = std::sqrt(x_norm);
    y_norm = std::sqrt(y_norm);
    z_norm = std::sqrt(z_norm);
    if (x == x_norm || y == y_norm || z == z_norm)
        return 2;
    else
        return 3;
}

void mdl_sld::write_prj_file(std::string& name) {
    std::ofstream sld_out_file(std::string(name + ".fes").c_str(),
                               std::ios::out | std::ios::ate | std::ios::app);
    sld_out_file << "#Sld_Nodes " << nodes.size() << "\n";
    for (size_t i = 0; i < nodes.size(); i++) {
        sld_out_file << std::scientific << std::setprecision(16) << nodes[i][0]
                     << " " << std::setprecision(16) << nodes[i][1] << " "
                     << std::setprecision(16) << nodes[i][2] << "\n";
    }
    if (dim == 2) {
        sld_out_file << "#Sld_Edges " << edges.size() << " "
                     << edges_marker[0].size() << "\n";
        for (size_t i = 0; i < edges.size(); i++) {
            for (size_t j = 0; j < edges[i].size(); j++)
                sld_out_file << edges[i][j] << " ";
            for (size_t j = 0; j < edges_marker[i].size(); j++)
                sld_out_file << edges_marker[i][j] << " ";
            sld_out_file << "\n";
        }
    }
    if (dim == 3) {
        sld_out_file << "#Sld_Faces " << faces.size() << "\n";
        for (size_t i = 0; i < faces.size(); i++) {
            sld_out_file << faces[i].polygons.size() << " "
                         << faces[i].holes.size() << " " << faces_marker[i].back()
                         << "\n";
            for (size_t j = 0; j < faces[i].polygons.size(); j++) {
                sld_out_file << faces[i].polygons[j].size() << " ";
                for (size_t k = 0; k < faces[i].polygons[j].size(); k++) {
                    sld_out_file << faces[i].polygons[j][k] << " ";
                }
                sld_out_file << "\n";
            }
            for (size_t j = 0; j < faces[i].holes.size(); j++) {
                sld_out_file << std::scientific << std::setprecision(16)
                             << faces[i].holes[j][0] << " " << std::setprecision(16)
                             << faces[i].holes[j][1] << " " << std::setprecision(16)
                             << faces[i].holes[j][2] << "\n";
            }
        }
    }
    sld_out_file << "#Sld_Holes " << holes.size() << "\n";
    for (size_t i = 0; i < holes.size(); i++) {
        sld_out_file << std::scientific << std::setprecision(16) << holes[i][0]
                     << " " << std::setprecision(16) << holes[i][1] << " "
                     << std::setprecision(16) << holes[i][2] << "\n";
    }
    sld_out_file << "#Sld_Regions " << regions.size() << "\n";
    for (size_t i = 0; i < regions.size(); i++) {
        sld_out_file << std::scientific << std::setprecision(16)
                     << regions[i][0] << " " << std::setprecision(16)
                     << regions[i][1] << " " << std::setprecision(16)
                     << regions[i][2] << " " << regions_marker[i] << "\n";
    }
    sld_out_file.close();
}

void mdl_sld::read_prj_file(std::string& name) {
    clear();
    std::ifstream sld_in_file(std::string(name + ".fes").c_str(), std::ios::in);
    std::string line;
    std::istringstream iss;
    unsigned int tmp_uint;
    double tmp_dbl;
    std::string tmp_str;
    if (sld_in_file.is_open()) {
        while (getline(sld_in_file, line)) {
            iss.clear();
            iss.str(line);
            iss >> tmp_str;
            if (strcmp(tmp_str.data(), "#Sld_Nodes") == 0) {
                iss >> tmp_uint;
                nodes.resize(tmp_uint);
                for (size_t i = 0; i < nodes.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    nodes[i].resize(3);
                    iss >> nodes[i][0];
                    iss >> nodes[i][1];
                    iss >> nodes[i][2];
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Edges") == 0) {
                dim = 2;
                iss >> tmp_uint;
                edges.resize(tmp_uint);
                edges_marker.resize(tmp_uint);
                iss >> tmp_uint;
                for (size_t i = 0; i < edges.size(); i++) {
                    edges[i].resize(2);
                    edges_marker[i].resize(tmp_uint);
                }
                for (size_t i = 0; i < edges.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> edges[i][0];
                    iss >> edges[i][1];
                    for (unsigned int j = 0; j < edges_marker[i].size(); j++) {
                        iss >> edges_marker[i][j];
                        max_edges_marker = std::max(max_edges_marker,
                                                    edges_marker[i][j]);
                        if (std::find(bc_markers.begin(), bc_markers.end(),
                                      edges_marker[i][j]) == bc_markers.end())
                            bc_markers.push_back(edges_marker[i][j]);
                    }
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Faces") == 0) {
                dim = 3;
                iss >> tmp_uint;
                faces.resize(tmp_uint);
                faces_marker.resize(tmp_uint);
                for (size_t i = 0; i < faces.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    unsigned int polygons, hole;
                    int bmark;
                    iss >> polygons;
                    iss >> hole;
                    iss >> bmark;
                    faces[i].polygons.resize(polygons);
                    faces[i].holes.resize(hole);
                    faces_marker[i].push_back(bmark);
                    max_faces_marker = std::max(max_faces_marker, bmark);
                    if (std::find(bc_markers.begin(), bc_markers.end(), bmark)
                            == bc_markers.end())
                        bc_markers.push_back(bmark);
                    for (unsigned int j = 0; j < faces[i].polygons.size();
                            j++) {
                        getline(sld_in_file, line);
                        iss.clear();
                        iss.str(line);
                        iss >> tmp_uint;
                        faces[i].polygons[j].resize(tmp_uint);
                        for (unsigned int k = 0;
                                k < faces[i].polygons[j].size(); k++)
                            iss >> faces[i].polygons[j][k];
                    }
                    for (unsigned int j = 0; j < faces[i].holes.size(); j++) {
                        getline(sld_in_file, line);
                        iss.clear();
                        iss.str(line);
                        faces[i].holes[j].resize(3);
                        for (unsigned int k = 0; k < 3; k++)
                            iss >> faces[i].holes[j][k];
                    }
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Holes") == 0) {
                iss >> tmp_uint;
                holes.resize(tmp_uint);
                for (size_t i = 0; i < holes.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    nodes[i].resize(3);
                    iss >> nodes[i][0];
                    iss >> nodes[i][1];
                    iss >> nodes[i][2];
                }
            }
            if (strcmp(tmp_str.data(), "#Sld_Regions") == 0) {
                iss >> tmp_uint;
                regions.resize(tmp_uint);
                regions_marker.resize(tmp_uint);
                for (size_t i = 0; i < regions.size(); i++) {
                    getline(sld_in_file, line);
                    iss.clear();
                    iss.str(line);
                    regions[i].resize(3);
                    iss >> regions[i][0];
                    iss >> regions[i][1];
                    iss >> regions[i][2];
                    iss >> regions_marker[i];
                    max_regions_marker = std::max(max_regions_marker,
                                                  regions_marker[i]);
                    if (std::find(mtrl_markers.begin(), mtrl_markers.end(),
                                  regions_marker[i]) == mtrl_markers.end())
                        mtrl_markers.push_back(regions_marker[i]);
                }
            }
        }
    }
    sld_in_file.close();
}

mdl_bc::mdl_bc() {
    label = 0;
}

mdl_bc::~mdl_bc() {
    std::vector<std::complex<double> >().swap(mode_beta);
    std::vector<std::vector<std::complex<double> > >().swap(mode_eig_vec);
    std::vector<std::vector<std::complex<double> > >().swap(mode_eig_vec_f);
    std::vector<std::vector<size_t> >().swap(mode_dof_map);
}

mdl_mtrl::mdl_mtrl() {
}

mdl_mtrl::~mdl_mtrl() {
}

void mdl_mtrl::upd_mtrl() {
    epsr2 = -tand * epsr;
}

void mdl_mtrl::upd_mtrl(double& freq) {
    epsr2 = -sigma / (2.0 * phys_const::pi * freq * phys_const::eps0)
            - tand * epsr;
}

double mdl_mtrl::calc_epsr2(double& freq) {
    return (-sigma / (2.0 * phys_const::pi * freq * phys_const::eps0)
            - tand * epsr);
}

mdl_msh::mdl_msh() {
    n_domains = 0;
    n_edges = 0;
    n_faces = 0;
    n_nodes = 0;
    n_tetras = 0;
}

void mdl_msh::write_prj_file(std::string& name) {
    std::ofstream msh_out_file(std::string(name + ".fes").c_str(),
                               std::ios::out | std::ios::ate | std::ios::app);
    msh_out_file << "#Mesh " << type << "\n";
    msh_out_file << "#Nodes " << n_nodes << "\n";
    for (size_t i = 0; i < n_nodes; i++) {
        msh_out_file << std::scientific << std::setprecision(16)
                     << nod_pos[i][0] << " " << std::setprecision(16)
                     << nod_pos[i][1] << " " << std::setprecision(16)
                     << nod_pos[i][2] << "\n";
    }
    msh_out_file << "#Edges " << n_edges << "\n";
    for (size_t i = 0; i < n_edges; i++) {
        msh_out_file << edg_nodes[i][0] << " " << edg_nodes[i][1] << " "
                     << edg_lab[i] << "\n";
    }
    msh_out_file << "#Faces " << n_faces << "\n";
    for (size_t i = 0; i < n_faces; i++) {
        msh_out_file << fac_nodes[i][0] << " " << fac_nodes[i][1] << " "
                     << fac_nodes[i][2] << " " << fac_lab[i] << "\n";
    }
    msh_out_file << "#Tetras " << n_tetras << "\n";
    for (size_t i = 0; i < n_tetras; i++) {
        msh_out_file << tet_nodes[i][0] << " " << tet_nodes[i][1] << " "
                     << tet_nodes[i][2] << " " << tet_nodes[i][3] << " "
                     << tet_lab[i] << "\n";
    }
    msh_out_file.close();
}

void mdl_msh::read_prj_file(std::string& name) {
    clear();
    std::ifstream msh_in_file(std::string(name + ".fes").c_str(), std::ios::in);
    std::string line;
    std::istringstream iss;
    unsigned int tmp_uint;
    double tmp_dbl;
    std::string tmp_str;
    if (msh_in_file.is_open()) {
        while (getline(msh_in_file, line)) {
            iss.clear();
            iss.str(line);
            iss >> tmp_str;
            if (strcmp(tmp_str.data(), "#Mesh") == 0) {
                iss >> type;
            }
            if (strcmp(tmp_str.data(), "#Nodes") == 0) {
                iss >> n_nodes;
                nod_pos.resize(n_nodes);
                for (size_t i = 0; i < n_nodes; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    nod_pos[i].resize(3);
                    iss >> nod_pos[i][0];
                    iss >> nod_pos[i][1];
                    iss >> nod_pos[i][2];
                }
            }
            if (strcmp(tmp_str.data(), "#Edges") == 0) {
                iss >> n_edges;
                edg_nodes.resize(n_edges);
                edg_lab.resize(n_edges);
                for (size_t i = 0; i < n_edges; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    edg_nodes[i].resize(2);
                    iss >> edg_nodes[i][0];
                    iss >> edg_nodes[i][1];
                    iss >> edg_lab[i];
                }
            }
            if (strcmp(tmp_str.data(), "#Faces") == 0) {
                iss >> n_faces;
                fac_nodes.resize(n_faces);
                fac_lab.resize(n_faces);
                for (size_t i = 0; i < n_faces; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    fac_nodes[i].resize(3);
                    iss >> fac_nodes[i][0];
                    iss >> fac_nodes[i][1];
                    iss >> fac_nodes[i][2];
                    iss >> fac_lab[i];
                }
            }
            if (strcmp(tmp_str.data(), "#Tetras") == 0) {
                iss >> n_tetras;
                tet_nodes.resize(n_tetras);
                tet_lab.resize(n_tetras);
                for (size_t i = 0; i < n_tetras; i++) {
                    getline(msh_in_file, line);
                    iss.clear();
                    iss.str(line);
                    tet_nodes[i].resize(4);
                    iss >> tet_nodes[i][0];
                    iss >> tet_nodes[i][1];
                    iss >> tet_nodes[i][2];
                    iss >> tet_nodes[i][3];
                    iss >> tet_lab[i];
                }
            }
        }
    }
    msh_in_file.close();
    regularize_mesh();
}

void mdl_msh::read_tetgen_files(std::string& name) {
    clear();
    type = "TETRA";
    unsigned int lvl = 1;
    std::ostringstream str_lvl;
    str_lvl << lvl;
    std::string line;
    std::istringstream iss;
    double tmp_dbl;
    int tmp_int;
    size_t tmp_uint;
    std::string tmp_str;
    std::ifstream tet_node_file(
        std::string(name + "." + str_lvl.str() + ".node").c_str());
    if (tet_node_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".node") << "\n";
        iss.clear();
        do {
            getline(tet_node_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_nodes;
        iss >> tmp_int;
        nod_pos.resize(n_nodes);
        for (size_t i = 0; i < n_nodes; i++) {
            iss.clear();
            do {
                getline(tet_node_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_dbl;
                nod_pos[i].push_back(tmp_dbl);
            }
        }
    }
    tet_node_file.close();
    std::ifstream tet_edge_file(
        std::string(name + "." + str_lvl.str() + ".edge").c_str());
    if (tet_edge_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".edge") << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tet_edge_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_edges;
        iss >> tmp_int;
        edg_nodes.resize(n_edges);
        edg_lab.assign(n_edges, 0);
        for (size_t i = 0; i < n_edges; i++) {
            iss.clear();
            do {
                getline(tet_edge_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 2; j++) {
                iss >> tmp_uint;
                edg_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            edg_lab[i] = tmp_int;
        }
    }
    tet_edge_file.close();
    std::ifstream tet_face_file(
        std::string(name + "." + str_lvl.str() + ".face").c_str());
    if (tet_face_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".face") << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tet_face_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_faces;
        iss >> tmp_int;
        fac_nodes.resize(n_faces);
        fac_lab.assign(n_faces, 0);
        for (size_t i = 0; i < n_faces; i++) {
            iss.clear();
            do {
                getline(tet_face_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_uint;
                fac_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> fac_lab[i];
        }
    }
    tet_face_file.close();
    std::ifstream tet_ele_file(
        std::string(name + "." + str_lvl.str() + ".ele").c_str());
    if (tet_ele_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".ele") << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tet_ele_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_tetras;
        iss >> tmp_int;
        tet_nodes.resize(n_tetras);
        tet_lab.assign(n_tetras, 0);
        for (size_t i = 0; i < n_tetras; i++) {
            iss.clear();
            do {
                getline(tet_ele_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 4; j++) {
                iss >> tmp_uint;
                tet_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            tet_lab[i] = tmp_int;
        }
    }
    tet_ele_file.close();
    regularize_mesh();
}

void mdl_msh::read_triangle_files(std::string& name) {
    clear();
    type = "TRIA";
    unsigned int lvl = 1;
    std::ostringstream str_lvl;
    str_lvl << lvl;
    std::string line;
    std::istringstream iss;
    double tmp_dbl;
    int tmp_int;
    size_t tmp_uint;
    std::string tmp_str;
    std::ifstream tria_node_file(
        std::string(name + "." + str_lvl.str() + ".node").c_str());
    if (tria_node_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".node") << "\n";
        iss.clear();
        do {
            getline(tria_node_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_nodes;
        iss >> tmp_int;
        nod_pos.resize(n_nodes);
        for (size_t i = 0; i < n_nodes; i++) {
            iss.clear();
            do {
                getline(tria_node_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_dbl;
                nod_pos[i].push_back(tmp_dbl);
            }
        }
    }
    tria_node_file.close();
    std::ifstream tria_edge_file(
        std::string(name + "." + str_lvl.str() + ".edge").c_str());
    if (tria_edge_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".edge") << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tria_edge_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_edges;
        iss >> tmp_int;
        edg_nodes.resize(n_edges);
        edg_lab.assign(n_edges, 0);
        for (size_t i = 0; i < n_edges; i++) {
            iss.clear();
            do {
                getline(tria_edge_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 2; j++) {
                iss >> tmp_uint;
                edg_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            edg_lab[i] = tmp_int;
        }
    }
    tria_edge_file.close();
    std::ifstream tria_ele_file(
        std::string(name + "." + str_lvl.str() + ".ele").c_str());
    if (tria_ele_file.is_open()) {
        std::cout << "Loading "
                  << std::string(name + "." + str_lvl.str() + ".ele") << "\n";
        unsigned int dim;
        iss.clear();
        do {
            getline(tria_ele_file, line);
        } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
        iss.str(line);
        iss >> n_faces;
        iss >> tmp_int;
        fac_nodes.resize(n_faces);
        fac_lab.assign(n_faces, 0);
        for (size_t i = 0; i < n_faces; i++) {
            iss.clear();
            do {
                getline(tria_ele_file, line);
            } while (line.compare(0, 1, "#") == 0 || line.size() == 0);
            iss.str(line);
            iss >> tmp_int;
            for (unsigned int j = 0; j < 3; j++) {
                iss >> tmp_uint;
                fac_nodes[i].push_back(tmp_uint - 1);
            }
            iss >> tmp_int;
            fac_lab[i] = tmp_int;
        }
    }
    tria_ele_file.close();
    regularize_mesh();
}

/// Refinement mapping
static unsigned int tet_tet_nodes_loc_map_id[] = { 0, 4, 5, 6, 1, 4, 7, 8, 2, 5,
                                                   7, 9, 3, 6, 8, 9, 4, 5, 6, 7, 4, 6, 7, 8, 5, 6, 7, 9, 6, 7, 8, 9
                                                 };

static unsigned int tet_fac_nodes_loc_map_id[] = { 1, 7, 8, 2, 7, 9, 3, 8, 9, 7,
                                                   8, 9, 0, 5, 6, 2, 5, 9, 3, 6, 9, 5, 6, 9, 0, 4, 6, 1, 4, 8, 3, 6, 8, 4,
                                                   6, 8, 0, 4, 5, 1, 4, 7, 2, 5, 7, 4, 5, 7, 4, 5, 6, 4, 6, 7, 4, 7, 8, 5,
                                                   6, 7, 5, 7, 9, 6, 7, 8, 6, 7, 9, 6, 8, 9
                                                 };

static unsigned int fac_fac_nodes_loc_map_id[] = { 0, 4, 5, 1, 3, 5, 2, 3, 4, 3,
                                                   4, 5
                                                 };

static unsigned int tet_edg_nodes_loc_map_id[] = { 0, 4, 0, 5, 0, 6, 1, 4, 1, 7,
                                                   1, 8, 2, 5, 2, 7, 2, 9, 3, 6, 3, 8, 3, 9, 4, 5, 4, 6, 4, 7, 4, 8, 5, 6,
                                                   5, 7, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9
                                                 };

static unsigned int fac_edg_nodes_loc_map_id[] = { 1, 3, 2, 3, 0, 4, 2, 4, 0, 5,
                                                   1, 5, 3, 4, 3, 5, 4, 5
                                                 };

static int tet_fac_lab_loc_map_id[] = { 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3,
                                        3, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1
                                      };

static int fac_edg_lab_loc_map_id[] = { 0, 0, 1, 1, 2, 2, -1, -1, -1 };

mdl_msh::mdl_msh(mdl_msh * msh) {
    msh->tet_nodes = tet_nodes;
    msh->tet_edges = tet_edges;
    msh->tet_faces = tet_faces;
    msh->fac_nodes = fac_nodes;
    msh->fac_edges = fac_edges;
    msh->edg_nodes = edg_nodes;
    msh->nod_pos = nod_pos;
    msh->edg_lab = edg_lab;
    msh->fac_lab = fac_lab;
    msh->tet_lab = tet_lab;
    msh->fac_adj_tet = fac_adj_tet;
    msh->edg_adj_fac = edg_adj_fac;
    msh->dom_tetras = dom_tetras;
    msh->dom_faces = dom_faces;
    msh->n_nodes = n_nodes;
    msh->n_edges = n_edges;
    msh->n_faces = n_faces;
    msh->n_tetras = n_tetras;
    msh->n_domains = n_domains;
}

mdl_msh::~mdl_msh() {
    clear();
}

void mdl_msh::get_mesh_statistics() {
    std::cout << "tet_nodes = " << tet_nodes.size() << "\n";
    std::cout << "tet_edges = " << tet_edges.size() << "\n";
    std::cout << "tet_faces = " << tet_faces.size() << "\n";
    std::cout << "fac_nodes = " << fac_nodes.size() << "\n";
    std::cout << "fac_edges = " << fac_edges.size() << "\n";
    std::cout << "edg_nodes = " << edg_nodes.size() << "\n";
    std::cout << "nod_pos = " << nod_pos.size() << "\n";
    std::cout << "fac_lab = " << fac_lab.size() << "\n";
    std::cout << "tet_lab = " << tet_lab.size() << "\n";
    std::cout << "fac_adj_tet = " << fac_adj_tet.size() << "\n";
    std::cout << "edg_adj_fac = " << edg_adj_fac.size() << "\n";
    std::cout << "dom_tetras = " << dom_tetras.size() << "\n";
    std::cout << "dom_faces = " << dom_faces.size() << "\n";
    std::cout << "n_nodes = " << n_nodes << "\n";
    std::cout << "n_edges = " << n_edges << "\n";
    std::cout << "n_faces = " << n_faces << "\n";
    std::cout << "n_tetras = " << n_tetras << "\n";
    std::cout << "n_domains = " << n_domains << "\n";
}

void mdl_msh::clear() {
    tet_nodes.clear();
    tet_edges.clear();
    tet_faces.clear();
    fac_nodes.clear();
    fac_edges.clear();
    edg_nodes.clear();
    nod_pos.clear();
    fac_lab.clear();
    tet_lab.clear();
    fac_adj_tet.clear();
    edg_adj_fac.clear();
    dom_tetras.clear();
    dom_faces.clear();
    n_nodes = 0;
    n_edges = 0;
    n_faces = 0;
    n_tetras = 0;
    n_domains = 0;
    max_edg_marker = -INT_MAX;
    max_fac_marker = -INT_MAX;
    max_tet_marker = -INT_MAX;
}

void mdl_msh::regularize_mesh() {
    max_edg_marker = -INT_MAX;
    max_fac_marker = -INT_MAX;
    max_tet_marker = -INT_MAX;
    std::map<std::pair<size_t, size_t>, size_t> edgesMap;
    std::map<std::tuple<size_t, size_t, size_t>, size_t> facesMap;
    if (strcmp(type.data(), "EDGE")) { // if not edges
        for (size_t i = 0; i < n_edges; i++) {
            std::sort(edg_nodes[i].begin(), edg_nodes[i].end());
            edgesMap[std::make_pair(edg_nodes[i][0], edg_nodes[i][1])] = i;
        }
        fac_edges.clear();
        edg_adj_fac.clear();
        fac_edges.resize(n_faces);
        edg_adj_fac.resize(n_edges);
        for (size_t i = 0; i < n_faces; i++) {
            std::sort(fac_nodes[i].begin(), fac_nodes[i].end());
            fac_edges[i].push_back(
                edgesMap[std::make_pair(fac_nodes[i][1], fac_nodes[i][2])]);
            fac_edges[i].push_back(
                edgesMap[std::make_pair(fac_nodes[i][0], fac_nodes[i][2])]);
            fac_edges[i].push_back(
                edgesMap[std::make_pair(fac_nodes[i][0], fac_nodes[i][1])]);
            facesMap[std::make_tuple(fac_nodes[i][0], fac_nodes[i][1],
                                                                      fac_nodes[i][2])] = i;
            for (size_t j = 0; j < 3; j++) {
                edg_adj_fac[fac_edges[i][j]].push_back(i);
            }
        }
    }
    if (strcmp(type.data(), "TETRA") == 0) {
        tet_edges.clear();
        tet_faces.clear();
        fac_adj_tet.clear();
        tet_edges.resize(n_tetras);
        tet_faces.resize(n_tetras);
        fac_adj_tet.resize(n_faces);
        for (size_t i = 0; i < n_tetras; i++) {
            std::sort(tet_nodes[i].begin(), tet_nodes[i].end());
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][0], tet_nodes[i][1])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][0], tet_nodes[i][2])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][0], tet_nodes[i][3])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][1], tet_nodes[i][2])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][1], tet_nodes[i][3])]);
            tet_edges[i].push_back(
                edgesMap[std::make_pair(tet_nodes[i][2], tet_nodes[i][3])]);
            tet_faces[i].push_back(
                facesMap[std::make_tuple(tet_nodes[i][1], tet_nodes[i][2],
                                                                          tet_nodes[i][3])]);
            tet_faces[i].push_back(
                facesMap[std::make_tuple(tet_nodes[i][0], tet_nodes[i][2],
                                                                          tet_nodes[i][3])]);
            tet_faces[i].push_back(
                facesMap[std::make_tuple(tet_nodes[i][0], tet_nodes[i][1],
                                                                          tet_nodes[i][3])]);
            tet_faces[i].push_back(
                facesMap[std::make_tuple(tet_nodes[i][0], tet_nodes[i][1],
                                                                          tet_nodes[i][2])]);
            for (size_t j = 0; j < 4; j++) {
                fac_adj_tet[tet_faces[i][j]].push_back(i);
            }
        }
    }
    for (size_t i = 0; i < edg_lab.size(); i++)
        max_edg_marker = std::max(max_edg_marker, edg_lab[i]);
    for (size_t i = 0; i < fac_lab.size(); i++)
        max_fac_marker = std::max(max_fac_marker, fac_lab[i]);
    for (size_t i = 0; i < tet_lab.size(); i++)
        max_tet_marker = std::max(max_tet_marker, tet_lab[i]);
}

void mdl_msh::refine_homogeneous() {
    std::vector<std::vector<double> > new_nod_pos(n_nodes + n_edges,
            std::vector<double>(3));
    std::vector<std::vector<size_t> > new_tet_nodes(n_tetras * 8,
            std::vector<size_t>(4));
    std::vector<std::vector<size_t> > new_tet_edges(n_tetras * 8,
            std::vector<size_t>(6));
    std::vector<std::vector<size_t> > new_tet_faces(n_tetras * 8,
            std::vector<size_t>(4));
    std::vector<std::vector<size_t> > new_fac_nodes(n_faces * 4 + n_tetras * 8,
            std::vector<size_t>(3));
    std::vector<std::vector<size_t> > new_fac_edges(n_faces * 4 + n_tetras * 8,
            std::vector<size_t>(3));
    std::vector<std::vector<size_t> > new_edg_nodes(
        n_edges * 2 + n_faces * 3 + n_tetras, std::vector<size_t>(2));
    std::vector<int> new_edg_lab(n_edges * 2 + n_faces * 3 + n_tetras, 0);
    std::vector<int> new_fac_lab(n_faces * 4 + n_tetras * 8, 0);
    std::vector<int> new_tet_lab(n_tetras * 8, 0);
    std::vector<std::vector<size_t> > new_edg_adj_fac(
        n_edges * 2 + n_faces * 3 + n_tetras);
    std::vector<std::vector<size_t> > new_fac_adj_tet(
        n_faces * 4 + n_tetras * 8);
//
    std::map<std::pair<size_t, size_t>, size_t> edgesMap;
    std::map<std::tuple<size_t, size_t, size_t>, size_t> facesMap;
//
    for (size_t nid = 0; nid < n_nodes; nid++)
        new_nod_pos[nid] = nod_pos[nid];
    for (size_t eid = 0; eid < n_edges; eid++) {
        new_nod_pos[n_nodes + eid][0] = (nod_pos[edg_nodes[eid][0]][0]
                                         + nod_pos[edg_nodes[eid][1]][0]) / 2;
        new_nod_pos[n_nodes + eid][1] = (nod_pos[edg_nodes[eid][0]][1]
                                         + nod_pos[edg_nodes[eid][1]][1]) / 2;
        new_nod_pos[n_nodes + eid][2] = (nod_pos[edg_nodes[eid][0]][2]
                                         + nod_pos[edg_nodes[eid][1]][2]) / 2;
    }
    size_t tet_lvl = 0;
    size_t fac_lvl = 0;
    size_t edg_lvl = 0;
    if (strcmp(type.data(), "TRIA") == 0) {
        std::vector<size_t> nod_glob(6), edg_tmp(2), fac_tmp(3);
        for (size_t fid = 0; fid < n_faces; fid++) {
            for (unsigned int i = 0; i < 3; i++) {
                nod_glob[i] = fac_nodes[fid][i];
                nod_glob[3 + i] = n_nodes + fac_edges[fid][i];
            }
            // populating edges
            for (unsigned int i = 0; i < 9; i++) {
                edg_tmp[0] = nod_glob[fac_edg_nodes_loc_map_id[2 * i]];
                edg_tmp[1] = nod_glob[fac_edg_nodes_loc_map_id[2 * i + 1]];
                std::sort(edg_tmp.begin(), edg_tmp.end());
                if (edgesMap.find(std::make_pair(edg_tmp[0], edg_tmp[1]))
                        == edgesMap.end()) {
                    edgesMap[std::make_pair(edg_tmp[0], edg_tmp[1])] = edg_lvl;
                    new_edg_nodes[edg_lvl][0] = edg_tmp[0];
                    new_edg_nodes[edg_lvl][1] = edg_tmp[1];
                    if (fac_edg_lab_loc_map_id[i] > -1) {
                        new_edg_lab[edg_lvl] =
                            edg_lab[fac_edges[fid][fac_edg_lab_loc_map_id[i]]];
                    }
                    ++edg_lvl;
                }
            }
            // populating faces
            for (unsigned int i = 0; i < 4; i++) {
                fac_tmp[0] = nod_glob[fac_fac_nodes_loc_map_id[3 * i]];
                fac_tmp[1] = nod_glob[fac_fac_nodes_loc_map_id[3 * i + 1]];
                fac_tmp[2] = nod_glob[fac_fac_nodes_loc_map_id[3 * i + 2]];
                std::sort(fac_tmp.begin(), fac_tmp.end());
                new_fac_nodes[fac_lvl][0] = fac_tmp[0];
                new_fac_nodes[fac_lvl][1] = fac_tmp[1];
                new_fac_nodes[fac_lvl][2] = fac_tmp[2];
                new_fac_edges[fac_lvl][0] = edgesMap[std::make_pair(fac_tmp[1],
                                                     fac_tmp[2])];
                new_fac_edges[fac_lvl][1] = edgesMap[std::make_pair(fac_tmp[0],
                                                     fac_tmp[2])];
                new_fac_edges[fac_lvl][2] = edgesMap[std::make_pair(fac_tmp[0],
                                                     fac_tmp[1])];
                new_fac_lab[fac_lvl] = fac_lab[fid];
                for (size_t j = 0; j < 3; j++) {
                    size_t eid = new_fac_edges[fac_lvl][j];
                    new_edg_adj_fac[eid].push_back(fac_lvl);
                }
                ++fac_lvl;
            }
        }
    }
    if (strcmp(type.data(), "TETRA") == 0) {
        std::vector<size_t> nod_glob(10), edg_tmp(2), fac_tmp(3), tet_tmp(4);
        for (size_t tid = 0; tid < n_tetras; tid++) {
            for (unsigned int i = 0; i < 4; i++)
                nod_glob[i] = tet_nodes[tid][i];
            for (unsigned int i = 0; i < 6; i++)
                nod_glob[4 + i] = n_nodes + tet_edges[tid][i];
            // populating edges
            for (unsigned int i = 0; i < 25; i++) {
                edg_tmp[0] = nod_glob[tet_edg_nodes_loc_map_id[2 * i]];
                edg_tmp[1] = nod_glob[tet_edg_nodes_loc_map_id[2 * i + 1]];
                std::sort(edg_tmp.begin(), edg_tmp.end());
                if (edgesMap.find(std::make_pair(edg_tmp[0], edg_tmp[1]))
                        == edgesMap.end()) {
                    edgesMap[std::make_pair(edg_tmp[0], edg_tmp[1])] = edg_lvl;
                    new_edg_nodes[edg_lvl][0] = edg_tmp[0];
                    new_edg_nodes[edg_lvl][1] = edg_tmp[1];
                    ++edg_lvl;
                }
            }
            // populating faces
            for (unsigned int i = 0; i < 24; i++) {
                fac_tmp[0] = nod_glob[tet_fac_nodes_loc_map_id[3 * i]];
                fac_tmp[1] = nod_glob[tet_fac_nodes_loc_map_id[3 * i + 1]];
                fac_tmp[2] = nod_glob[tet_fac_nodes_loc_map_id[3 * i + 2]];
                std::sort(fac_tmp.begin(), fac_tmp.end());
                if (facesMap.find(
                            std::make_tuple(fac_tmp[0], fac_tmp[1], fac_tmp[2]))
                        == facesMap.end()) {
                    facesMap[std::make_tuple(fac_tmp[0], fac_tmp[1], fac_tmp[2])] =
                        fac_lvl;
                    new_fac_nodes[fac_lvl][0] = fac_tmp[0];
                    new_fac_nodes[fac_lvl][1] = fac_tmp[1];
                    new_fac_nodes[fac_lvl][2] = fac_tmp[2];
                    new_fac_edges[fac_lvl][0] = edgesMap[std::make_pair(
                            fac_tmp[1], fac_tmp[2])];
                    new_fac_edges[fac_lvl][1] = edgesMap[std::make_pair(
                            fac_tmp[0], fac_tmp[2])];
                    new_fac_edges[fac_lvl][2] = edgesMap[std::make_pair(
                            fac_tmp[0], fac_tmp[1])];
                    if (tet_fac_lab_loc_map_id[i] > -1) {
                        new_fac_lab[fac_lvl] =
                            fac_lab[tet_faces[tid][tet_fac_lab_loc_map_id[i]]];
                    }
                    ++fac_lvl;
                }
            }
            for (unsigned int i = 0; i < 8; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    tet_tmp[j] = nod_glob[tet_tet_nodes_loc_map_id[i * 4 + j]];
                }
                std::sort(tet_tmp.begin(), tet_tmp.end());
                for (unsigned int j = 0; j < 4; j++) {
                    new_tet_nodes[tet_lvl][j] = tet_tmp[j];
                }
                new_tet_lab[tet_lvl] = tet_lab[tid];
                new_tet_edges[tet_lvl][0] = edgesMap[std::make_pair(tet_tmp[0],
                                                     tet_tmp[1])];
                new_tet_edges[tet_lvl][1] = edgesMap[std::make_pair(tet_tmp[0],
                                                     tet_tmp[2])];
                new_tet_edges[tet_lvl][2] = edgesMap[std::make_pair(tet_tmp[0],
                                                     tet_tmp[3])];
                new_tet_edges[tet_lvl][3] = edgesMap[std::make_pair(tet_tmp[1],
                                                     tet_tmp[2])];
                new_tet_edges[tet_lvl][4] = edgesMap[std::make_pair(tet_tmp[1],
                                                     tet_tmp[3])];
                new_tet_edges[tet_lvl][5] = edgesMap[std::make_pair(tet_tmp[2],
                                                     tet_tmp[3])];
                new_tet_faces[tet_lvl][0] = facesMap[std::make_tuple(tet_tmp[1],
                                                     tet_tmp[2], tet_tmp[3])];
                new_tet_faces[tet_lvl][1] = facesMap[std::make_tuple(tet_tmp[0],
                                                     tet_tmp[2], tet_tmp[3])];
                new_tet_faces[tet_lvl][2] = facesMap[std::make_tuple(tet_tmp[0],
                                                     tet_tmp[1], tet_tmp[3])];
                new_tet_faces[tet_lvl][3] = facesMap[std::make_tuple(tet_tmp[0],
                                                     tet_tmp[1], tet_tmp[2])];
                for (size_t j = 0; j < 4; j++) {
                    size_t fid = new_tet_faces[tet_lvl][j];
                    new_fac_adj_tet[fid].push_back(tet_lvl);
                }
                ++tet_lvl;
            }
        }
    }
    nod_pos = new_nod_pos;
    edg_nodes = new_edg_nodes;
    edg_lab = new_edg_lab;
    edg_adj_fac = new_edg_adj_fac;
    fac_nodes = new_fac_nodes;
    fac_edges = new_fac_edges;
    fac_lab = new_fac_lab;
    fac_adj_tet = new_fac_adj_tet;
    tet_nodes = new_tet_nodes;
    tet_edges = new_tet_edges;
    tet_faces = new_tet_faces;
    tet_lab = new_tet_lab;
    n_nodes = n_nodes + n_edges;
    n_edges = n_edges * 2 + n_faces * 3 + n_tetras;
    n_faces = n_faces * 4 + n_tetras * 8;
    n_tetras = n_tetras * 8;
}

void mdl_msh::save_vtk_mesh(std::string vtkMshName) {
    unsigned int n_bc_mtrl_faces = 0;
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] > -1)
            ++n_bc_mtrl_faces;
        else if (fac_adj_tet[i].size() > 1) {
            if (tet_lab[fac_adj_tet[i][0]] != tet_lab[fac_adj_tet[i][1]]) {
                --fac_lab[i];
                ++n_bc_mtrl_faces;
            }
        }
    }
    std::ofstream out_vol_msh(std::string(vtkMshName + "_volmsh.vtk").data());
    out_vol_msh << "# vtk DataFile Version 2.0\n";
    out_vol_msh << "Mesh data\n";
    out_vol_msh << "ASCII\n";
    out_vol_msh << "DATASET UNSTRUCTURED_GRID\n";
    out_vol_msh << "POINTS " << n_nodes << " double \n";
    for (size_t i = 0; i < n_nodes; i++) {
        out_vol_msh << std::setprecision(16) << nod_pos[i][0] << " ";
        out_vol_msh << std::setprecision(16) << nod_pos[i][1] << " ";
        out_vol_msh << std::setprecision(16) << nod_pos[i][2] << "\n";
    }
    out_vol_msh << "CELLS " << n_tetras << " " << 5 * n_tetras << "\n";
    for (size_t i = 0; i < n_tetras; i++) {
        out_vol_msh << 4 << " ";
        out_vol_msh << tet_nodes[i][0] << " ";
        out_vol_msh << tet_nodes[i][1] << " ";
        out_vol_msh << tet_nodes[i][2] << " ";
        out_vol_msh << tet_nodes[i][3] << "\n";
    }
    out_vol_msh << "CELL_TYPES " << n_tetras << "\n";
    for (size_t i = 0; i < n_tetras; i++) {
        out_vol_msh << 10 << "\n";
    }
    out_vol_msh << "CELL_DATA " << n_tetras << "\n";
    out_vol_msh << "SCALARS " << "Materials int 1\n";
    out_vol_msh << "LOOKUP_TABLE jet\n";
    for (size_t i = 0; i < n_tetras; i++) {
        out_vol_msh << tet_lab[i] << "\n";
    }
    out_vol_msh.close();
    std::ofstream out_srf_msh(std::string(vtkMshName + "_srfmsh.vtk").data());
    out_srf_msh << "# vtk DataFile Version 2.0\n";
    out_srf_msh << "Mesh data\n";
    out_srf_msh << "ASCII\n";
    out_srf_msh << "DATASET UNSTRUCTURED_GRID\n";
    out_srf_msh << "POINTS " << n_nodes << " double \n";
    for (size_t i = 0; i < n_nodes; i++) {
        out_srf_msh << std::setprecision(16) << nod_pos[i][0] << " ";
        out_srf_msh << std::setprecision(16) << nod_pos[i][1] << " ";
        out_srf_msh << std::setprecision(16) << nod_pos[i][2] << "\n";
    }
    out_srf_msh << "CELLS " << n_bc_mtrl_faces << " " << 4 * n_bc_mtrl_faces
                << "\n";
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] != -1) {
            out_srf_msh << 3 << " ";
            out_srf_msh << fac_nodes[i][0] << " ";
            out_srf_msh << fac_nodes[i][1] << " ";
            out_srf_msh << fac_nodes[i][2] << "\n";
        }
    }
    out_srf_msh << "CELL_TYPES " << n_bc_mtrl_faces << "\n";
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] != -1)
            out_srf_msh << 5 << "\n";
    }
    out_srf_msh << "CELL_DATA " << n_bc_mtrl_faces << "\n";
    out_srf_msh << "SCALARS " << "Boundaries int 1\n";
    out_srf_msh << "LOOKUP_TABLE jet\n";
    for (size_t i = 0; i < n_faces; i++) {
        if (fac_lab[i] != -1)
            out_srf_msh << fac_lab[i] << "\n";
    }
    out_srf_msh.close();
    for (size_t i = 0; i < n_faces; i++)
        if (fac_lab[i] < -1)
            fac_lab[i] = -1;
}

std::vector<std::vector<double> > mdl_msh::tet_geo(size_t id) {
    std::vector<std::vector<double> > cGeo(4, std::vector<double>(3));
    cGeo[0] = nod_pos[tet_nodes[id][0]];
    cGeo[1] = nod_pos[tet_nodes[id][1]];
    cGeo[2] = nod_pos[tet_nodes[id][2]];
    cGeo[3] = nod_pos[tet_nodes[id][3]];
    return cGeo;
}

std::vector<std::vector<double> > mdl_msh::fac_geo(size_t id) {
    std::vector<std::vector<double> > cGeo(3, std::vector<double>(3));
    cGeo[0] = nod_pos[fac_nodes[id][0]];
    cGeo[1] = nod_pos[fac_nodes[id][1]];
    cGeo[2] = nod_pos[fac_nodes[id][2]];
    return cGeo;
}

std::vector<std::vector<double> > mdl_msh::fac_geo2(size_t id) {
    std::vector<double> v0 = nod_pos[fac_nodes[id][0]];
    std::vector<double> v1 = nod_pos[fac_nodes[id][1]];
    std::vector<double> v2 = nod_pos[fac_nodes[id][2]];
    for (unsigned int i = 0; i < 3; i++) {
        v1[i] -= v0[i];
        v2[i] -= v0[i];
    }
    std::vector<double> u = v1;
    std::vector<double> n(3), v(3);
    n[0] = v1[1] * v2[2];
    n[1] = v1[2] * v2[0];
    n[2] = v1[0] * v2[1];
    double norm2_u = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    double norm2_n = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (unsigned int i = 0; i < 3; i++) {
        u[i] /= norm2_u;
        n[i] /= norm2_n;
    }
    v[0] = n[1] * u[2];
    v[1] = n[2] * u[0];
    v[2] = n[0] * u[1];
    std::vector<std::vector<double> > cGeo2(3, std::vector<double>(2, 0.0));
    cGeo2[1][0] = v1[0] * u[0] + v1[1] * u[1] + v1[2] * u[2];
    cGeo2[2][0] = v2[0] * u[0] + v2[1] * u[1] + v2[2] * u[2];
    cGeo2[2][1] = v2[0] * v[0] + v2[1] * v[1] + v2[2] * v[2];
    return cGeo2;
}

std::vector<std::vector<double> > mdl_msh::edg_geo(size_t id) {
    std::vector<std::vector<double> > cGeo(2, std::vector<double>(3));
    cGeo[0] = nod_pos[edg_nodes[id][0]];
    cGeo[1] = nod_pos[edg_nodes[id][1]];
    return cGeo;
}

std::vector<double> mdl_msh::int_node(size_t id) {
    std::vector<size_t> nfac = fac_nodes[id];
    std::vector<size_t> ntet = tet_nodes[fac_adj_tet[id][0]];
    size_t intid;
    for (unsigned int i = 0; i < 4; i++) {
        bool found = true;
        intid = ntet[i];
        for (unsigned int j = 0; j < 3; j++)
            if (ntet[i] == nfac[j]) {
                found = false;
            }
        if (found) {
            break;
        }
    }
    return nod_pos[intid];
}

std::vector<double> mdl_msh::int_node(size_t id, size_t& ref_face) {
    std::vector<size_t> nfac = fac_nodes[id];
    std::vector<size_t> ntet = tet_nodes[fac_adj_tet[id][ref_face]];
    size_t intid = 0;
    for (unsigned int i = 0; i < 4; i++) {
        bool found = true;
        intid = ntet[i];
        for (unsigned int j = 0; j < 3; j++)
            if (ntet[i] == nfac[j]) {
                found = false;
            }
        if (found) {
            ref_face = i;
            break;
        }
    }
    return nod_pos[intid];
}

//double mdl_msh::tet_vol(size_t id) {
////         1      [ax bx cx dx]T
////    V = --- det [ay by cy dy]
////         6      [az bz cz dz]
////                [ 1  1  1  1]
//    std::vector<std::vector<double> > matrix = tet_geo(id);
//    matrix.resize(4,4);
//    matrix.col(3).fill(1.0);
//    return std::abs(arma::det(matrix]]/6.0;
//}
//
//double mdl_msh::tet_mean_edg_length(size_t id) {
//    std::vector<std::vector<double> > matrix = tet_geo(id);
//    return (arma::norm(matrix[1)-matrix[0),2) +
//            arma::norm(matrix[2)-matrix[0),2) +
//            arma::norm(matrix[3)-matrix[0),2) +
//            arma::norm(matrix[2)-matrix[1),2) +
//            arma::norm(matrix[3)-matrix[1),2) +
//            arma::norm(matrix[3)-matrix[2),2]] / 6.0;
//}
//
//double mdl_msh::tet_max_edg_length(size_t id) {
//    std::vector<std::vector<double> > matrix = tet_geo(id);
//    double max_edge = 0.0;
//    max_edge = std::max(max_edge, arma::norm(matrix[1)-matrix[0),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[2)-matrix[0),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[3)-matrix[0),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[2)-matrix[1),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[3)-matrix[1),2]];
//    max_edge = std::max(max_edge, arma::norm(matrix[3)-matrix[2),2]];
//    return max_edge;
//}

void mdl_msh::get_bounding_info() {
    bounding_box.resize(2);
    bounding_box[0].assign(3, DBL_MAX);
    bounding_box[1].assign(3, -DBL_MAX);
    for (size_t i = 0; i < nod_pos.size(); i++) {
        bounding_box[0][0] = std::min(bounding_box[0][0], nod_pos[i][0]);
        bounding_box[0][1] = std::min(bounding_box[0][1], nod_pos[i][1]);
        bounding_box[0][2] = std::min(bounding_box[0][2], nod_pos[i][2]);
        bounding_box[1][0] = std::max(bounding_box[1][0], nod_pos[i][0]);
        bounding_box[1][1] = std::max(bounding_box[1][1], nod_pos[i][1]);
        bounding_box[1][2] = std::max(bounding_box[1][2], nod_pos[i][2]);
    }
    double x_dim, y_dim, z_dim;
    x_dim = std::abs(bounding_box[1][0] - bounding_box[0][0]);
    y_dim = std::abs(bounding_box[1][1] - bounding_box[0][1]);
    z_dim = std::abs(bounding_box[1][2] - bounding_box[0][2]);
    max_dimension = std::max(max_dimension, x_dim);
    max_dimension = std::max(max_dimension, y_dim);
    max_dimension = std::max(max_dimension, z_dim);
}

double mdl_msh::get_geom_dim() {
    double dim = 0;
    for (unsigned int i = 0; i < nod_pos.size(); i++) {
        dim = std::max(dim, std::abs(nod_pos[i][0]));
        dim = std::max(dim, std::abs(nod_pos[i][1]));
        dim = std::max(dim, std::abs(nod_pos[i][2]));
    }
    dim *= 1e-1;
    for (int i = -15; i < 15; i++) {
        if (dim < pow(10, double(i)))
            return pow(10, double(i));
    }
    return 0;
}

mdl_src::mdl_src() {
}

mdl_src::~mdl_src() {
}

mdl_frm::mdl_frm() {
}

mdl_frm::~mdl_frm() {
    std::vector<mdl_bc>().swap(bcs);
    std::vector<mdl_mtrl>().swap(mtrls);
}

void mdl_frm::clear() {
    bcs.clear();
    mtrls.clear();
}

void mdl_frm::update_msh_info(mdl_msh& msh) {
    if (strcmp(type.data(), "TETRA") == 0) {
        for (size_t i = 0; i < bcs.size(); i++) {
            bcs[i].faces.clear();
            size_t cLab = bcs[i].label;
            for (size_t fid = 0; fid < msh.n_faces; fid++) {
                if (cLab == msh.fac_lab[fid]) {
                    bcs[i].faces.push_back(fid);
                }
            }
        }
        std::vector<size_t> mtrl_tets(mtrls.size(), 0), LabMap(mtrls.size(), 0);
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            mtrl_tets[msh.tet_lab[tid]]++;
        }
        for (size_t i = 0; i < mtrls.size(); i++) {
            mtrls[i].tetras.resize(mtrl_tets[i]);
            mtrl_tets[i] = 0;
        }
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            size_t cLab = msh.tet_lab[tid];
            mtrls[cLab].tetras[mtrl_tets[cLab]++] = tid;
        }
    }
}

void mdl_frm::write_prj_file(std::string& name) {
    std::ofstream prj_out_file(std::string(name + ".fes").c_str(),
                               std::ios::out | std::ios::ate);
    prj_out_file << "#Formulation " << type << "\n";
    prj_out_file << "#Materials " << mtrls.size() << "\n";
    for (size_t i = 0; i < mtrls.size(); i++) {
        prj_out_file << mtrls[i].label << " " << mtrls[i].type << " "
                     << mtrls[i].epsr << " " << mtrls[i].mur << " " << mtrls[i].sigma << " "
                     << mtrls[i].tand << " " << mtrls[i].name << "\n";
    }
    prj_out_file << "#Boundaries " << bcs.size() << "\n";
    for (size_t i = 0; i < bcs.size(); i++) {
        prj_out_file << bcs[i].label << " " << bcs[i].name << " "  << bcs[i].type;
        if (strcmp(bcs[i].type.data(), "WavePort") == 0) {
            prj_out_file << " " << bcs[i].num_modes;
        } else if (strcmp(bcs[i].type.data(), "Impedance") == 0) {
            prj_out_file << " " << bcs[i].surf_impedance.real() << " " << bcs[i].surf_impedance.imag();
        } else if (strcmp(bcs[i].type.data(), "LumpedPort") == 0) {
            prj_out_file << " " << bcs[i].lumped_impedance.real() << " " << bcs[i].lumped_impedance.imag();
        } else if (strcmp(bcs[i].type.data(), "LumpedRLC") == 0) {
            prj_out_file << " " << bcs[i].R << " " << bcs[i].L << " " << bcs[i].C;
        } else if (strcmp(bcs[i].type.data(), "Voltage") == 0) {
            prj_out_file << " " << bcs[i].voltage;
        } else if (strcmp(bcs[i].type.data(), "Current") == 0) {
            prj_out_file << " " << bcs[i].current;
        }
        prj_out_file << "\n";
    }
    prj_out_file.close();
}

void mdl_frm::read_prj_file(std::string& name) {
    clear();
    std::ifstream frm_in_file(std::string(name + ".fes").c_str(), std::ios::in);
    std::string line;
    std::istringstream iss;
    unsigned int tmp_uint;
    double tmp_dbl;
    std::string tmp_str;
    if (frm_in_file.is_open()) {
        while (getline(frm_in_file, line)) {
            iss.clear();
            iss.str(line);
            iss >> tmp_str;
            if (strcmp(tmp_str.data(), "#Formulation") == 0) {
                iss >> type;
            }
            if (strcmp(tmp_str.data(), "#Materials") == 0) {
                iss >> tmp_uint;
                mtrls.resize(tmp_uint);
                for (size_t i = 0; i < mtrls.size(); i++) {
                    getline(frm_in_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> mtrls[i].label;
                    iss >> mtrls[i].type;
                    iss >> mtrls[i].epsr;
                    iss >> mtrls[i].mur;
                    iss >> mtrls[i].sigma;
                    iss >> mtrls[i].tand;
                    iss >> mtrls[i].name;
                }
            }
            if (strcmp(tmp_str.data(), "#Boundaries") == 0) {
                iss >> tmp_uint;
                bcs.resize(tmp_uint);
                for (size_t i = 0; i < bcs.size(); i++) {
                    getline(frm_in_file, line);
                    iss.clear();
                    iss.str(line);
                    iss >> bcs[i].label;
                    iss >> bcs[i].name;
                    iss >> bcs[i].type;
                    if (strcmp(bcs[i].type.data(), "WavePort") == 0) {
                        iss >> bcs[i].num_modes;
                    } else if (strcmp(bcs[i].type.data(), "Impedance") == 0) {
                        double real;
                        double imag;
                        iss >> real;
                        iss >> imag;
                        bcs[i].surf_impedance = std::complex<double>(real,
                                                imag);
                    } else if (strcmp(bcs[i].type.data(), "LumpedPort") == 0) {
                        double real;
                        double imag;
                        iss >> real;
                        iss >> imag;
                        bcs[i].lumped_impedance = std::complex<double>(real, imag);
                    } else if (strcmp(bcs[i].type.data(), "LumpedRLC") == 0) {
                        iss >> bcs[i].R;
                        iss >> bcs[i].L;
                        iss >> bcs[i].C;
                    } else if (strcmp(bcs[i].type.data(), "Voltage") == 0) {
                        iss >> bcs[i].voltage;
                    } else if (strcmp(bcs[i].type.data(), "Current") == 0) {
                        iss >> bcs[i].current;
                    }
                }
            }
        }
    }
    frm_in_file.close();
}

void mdl_core::wrap_hfss(std::string& name, std::string& aux_path) {
    struct hfss_part {
        std::string name;
        std::string material;
        bool solve_inside = false;
        size_t id;
    };
    struct hfss_bnd {
        std::string name;
        std::string type;
        std::vector<size_t> faces;
        std::vector<size_t> solids;
        std::vector<size_t> face_ids;
        int num_modes;
        double FullResistance = 0.0;
        double FullReactance = 0.0;
        double Resistance = 0.0;
        double Reactance = 0.0;
        double Inductance = 0.0;
        double Capacitance = 0.0;
    };
    struct hfss_mtrl {
        double permittivity = 1.0;
        double permeability = 1.0;
        double conductivity = 0.0;
        double dielectric_loss_tangent = 0.0;
        std::string name;
    };
    bool debug = false;
    std::map<std::string, hfss_mtrl> mtrls;
    std::vector<size_t> mtrl_tag;
    std::vector<size_t> hfss_id;
    std::vector<bool> tet_flag;
    std::vector<bool> fac_flag;
    std::vector<bool> node_flag;
    std::vector<std::vector<size_t> > fac_hfss_tag;
    std::vector<hfss_bnd> bnds;
    std::vector<hfss_part> parts;
    std::map<size_t, std::vector<size_t> > bnd_map;
    std::vector<std::vector<size_t> > adj_tetra;

// starting algorithm
    msh.type = "TETRA";
    frm.type = "EM_E_FD";
    std::cout << "- Parsing HFSS project files\n";
//  READ MAIN HFSS
    {
        std::string R, L, C, X;
        double tmp_dbl;
        std::string tmp_str;
        std::string partName;
        // bool solve_inside;
        std::string materialName;
        std::string boundaryName;
        std::string line;
        std::ifstream file_name(std::string(name + ".hfss").c_str());
        if (file_name.is_open()) {
            while (getline(file_name, line)) { //file_name.good())
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "$begin") {
                    iss >> tmp_str;
                    if (tmp_str == "\'Materials\'") {
                        while (getline(file_name, line)) {
                            std::istringstream iss(line);
                            iss >> tmp_str;
                            if (tmp_str == "$begin") {
                                hfss_mtrl hfssMaterial;
                                std::string tmpMtrl;
                                while (iss.good()) {
                                    iss >> tmp_str;
                                    tmpMtrl.append(tmp_str);
                                }
                                materialName = tmpMtrl.substr(1,
                                                              tmpMtrl.size() - 2);
                                while (getline(file_name, line)) {
                                    std::istringstream iss(line);
                                    iss >> tmp_str; // $begin
                                    if (tmp_str.substr(0, 12)
                                            == "permittivity") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        14,
                                                        tmp_str.size()
                                                        - 15)).data());
                                        hfssMaterial.permittivity = tmp_dbl;
                                    } else if (tmp_str.substr(0, 12)
                                               == "permeability") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        14,
                                                        tmp_str.size()
                                                        - 15)).data());
                                        hfssMaterial.permeability = tmp_dbl;
                                    } else if (tmp_str.substr(0, 12)
                                               == "conductivity") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        14,
                                                        tmp_str.size()
                                                        - 15)).data());
                                        hfssMaterial.conductivity = tmp_dbl;
                                    } else if (tmp_str.substr(0, 23)
                                               == "dielectric_loss_tangent") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        25,
                                                        tmp_str.size()
                                                        - 26)).data());
                                        hfssMaterial.dielectric_loss_tangent =
                                            tmp_dbl;
                                    } else if (tmp_str == "$end") {
                                        std::string tmpMtrl;
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            tmpMtrl.append(tmp_str);
                                        }
                                        tmpMtrl = tmpMtrl.substr(1,
                                                                 tmpMtrl.size() - 2);
                                        if (tmpMtrl == materialName) {
                                            hfssMaterial.name = materialName;
                                            mtrls[materialName] = hfssMaterial;
                                            break;
                                        }
                                    }
                                }
                            } else if (tmp_str == "$end") {
                                iss >> tmp_str;
                                if (tmp_str == "\'Materials\'") {
                                    break;
                                }
                            }
                        }
                    } else if (tmp_str == "\'ToplevelParts\'") {
                        while (getline(file_name, line)) {
                            std::istringstream iss(line);
                            iss >> tmp_str;
                            if (tmp_str == "$begin") {
                                iss >> tmp_str;
                                if (tmp_str == "\'GeometryPart\'") {
                                    hfss_part hfssPart;
                                    while (getline(file_name, line)) {
                                        std::istringstream iss(line);
                                        iss >> tmp_str;
                                        if (tmp_str == "$begin") {
                                            iss >> tmp_str;
                                            if (tmp_str == "\'Attributes\'") {
                                                while (getline(file_name, line)) {
                                                    std::istringstream iss(
                                                        line);
                                                    iss >> tmp_str;
                                                    if (tmp_str.substr(0, 4)
                                                            == "Name") {
                                                        std::string tmpString =
                                                            tmp_str;
                                                        while (iss.good()) {
                                                            iss >> tmp_str;
                                                            tmpString.append(
                                                                tmp_str);
                                                        }
                                                        hfssPart.name =
                                                            std::string(
                                                                tmpString.substr(
                                                                    6,
                                                                    tmpString.size()
                                                                    - 7));
                                                    } else if (tmp_str.substr(0,
                                                                              13)
                                                               == "MaterialValue") { // for hfss v13
                                                        std::string tmpString =
                                                            tmp_str;
                                                        while (iss.good()) {
                                                            iss >> tmp_str;
                                                            tmpString.append(
                                                                tmp_str);
                                                        }
                                                        hfssPart.material =
                                                            std::string(
                                                                tmpString.substr(
                                                                    16,
                                                                    tmpString.size()
                                                                    - 18).data());
                                                    } else if (tmp_str.substr(0,
                                                                              12)
                                                               == "MaterialName") { // for hfss v11
                                                        std::string tmpString =
                                                            tmp_str;
                                                        while (iss.good()) {
                                                            iss >> tmp_str;
                                                            tmpString.append(
                                                                tmp_str);
                                                        }
                                                        hfssPart.material =
                                                            std::string(
                                                                tmpString.substr(
                                                                    14,
                                                                    tmpString.size()
                                                                    - 15));
                                                    } else if (tmp_str.substr(0,
                                                                              11)
                                                               == "SolveInside") {
                                                        hfssPart.solve_inside =
                                                            std::string(
                                                                tmp_str.substr(
                                                                    12,
                                                                    tmp_str.size()
                                                                    - 12))
                                                            == "true";
                                                    } else if (tmp_str
                                                               == "$end") {
                                                        iss >> tmp_str;
                                                        if (tmp_str
                                                                == "\'Attributes\'") {
                                                            break;
                                                        }
                                                    }
                                                }
                                            } else if (tmp_str
                                                       == "\'Operation\'") {
                                                while (getline(file_name, line)) {
                                                    std::istringstream iss(
                                                        line);
                                                    iss >> tmp_str;
                                                    if (tmp_str.substr(0, 12)
                                                            == "ParentPartID") {
                                                        hfssPart.id =
                                                            atoi(
                                                                std::string(
                                                                    tmp_str.substr(
                                                                        13,
                                                                        tmp_str.size()
                                                                        - 13)).data());
                                                    } else if (tmp_str
                                                               == "$end") {
                                                        iss >> tmp_str;
                                                        if (tmp_str
                                                                == "\'Operation\'") {
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        } else if (tmp_str == "$end") {
                                            iss >> tmp_str;
                                            if (tmp_str == "\'GeometryPart\'") {
                                                parts.push_back(hfssPart);
                                                break;
                                            }
                                        }
                                    }
                                }
                            } else if (tmp_str == "$end") {
                                iss >> tmp_str;
                                if (tmp_str == "\'ToplevelParts\'") {
                                    break;
                                }
                            }
                        }
                    } else if (tmp_str == "\'Boundaries\'") {
                        while (getline(file_name, line)) {
                            std::istringstream iss(line);
                            iss >> tmp_str; // $begin
                            if (tmp_str == "$begin") {
                                iss >> tmp_str;
                                hfss_bnd hfssBoundary;
                                hfssBoundary.name = tmp_str.substr(1,
                                                                   tmp_str.size() - 2);
                                while (getline(file_name, line)) {
                                    std::istringstream iss(line);
                                    iss >> tmp_str; // $begin
                                    if (tmp_str.substr(0, 9) == "BoundType") {
                                        hfssBoundary.type = std::string(
                                                                tmp_str.substr(11,
                                                                               tmp_str.size() - 11));
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            hfssBoundary.type += std::string(
                                                                     tmp_str.substr(0,
                                                                                    tmp_str.size()
                                                                                    - 1));
                                        }
                                        if (hfssBoundary.type.substr(0, 9)
                                                == "Radiation") {
                                            hfssBoundary.type = "Radiation";
                                        } else if (hfssBoundary.type.substr(0,
                                                                            9) == "Impedance") {
                                            hfssBoundary.type = "Impedance";
                                        }
                                    } else if (tmp_str.substr(0, 8)
                                               == "NumModes") {
                                        hfssBoundary.num_modes =
                                            atoi(
                                                std::string(
                                                    tmp_str.substr(
                                                        9,
                                                        tmp_str.size()
                                                        - 1)).data());
                                    } else if (tmp_str.substr(0, 14)
                                               == "FullResistance") {
                                        hfssBoundary.FullResistance =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        16,
                                                        tmp_str.size()
                                                        - 20)).data());
                                    } else if (tmp_str.substr(0, 13)
                                               == "FullReactance") {
                                        hfssBoundary.FullReactance =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        15,
                                                        tmp_str.size()
                                                        - 19)).data());
                                    } else if (tmp_str.substr(0, 10)
                                               == "Resistance") {
                                        R = tmp_str.substr(11,
                                                           tmp_str.size() - 11).data();
                                        removeCharsFromString(R, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(R));
                                        removeCharsFromString(R, SI_chars);
                                        hfssBoundary.Resistance = atof(R.data())
                                                                  * factor;
                                    } else if (tmp_str.substr(0, 9)
                                               == "Reactance") {
                                        X = tmp_str.substr(10,
                                                           tmp_str.size() - 10).data();
                                        removeCharsFromString(X, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(X));
                                        removeCharsFromString(X, SI_chars);
                                        hfssBoundary.Reactance = atof(X.data())
                                                                 * factor;
                                    } else if (tmp_str.substr(0, 10)
                                               == "Inductance") {
                                        L = tmp_str.substr(11,
                                                           tmp_str.size() - 11).data();
                                        removeCharsFromString(L, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(L));
                                        removeCharsFromString(L, SI_chars);
                                        hfssBoundary.Inductance = atof(L.data())
                                                                  * factor;
                                    } else if (tmp_str.substr(0, 11)
                                               == "Capacitance") {
                                        C = tmp_str.substr(12,
                                                           tmp_str.size() - 12).data();
                                        removeCharsFromString(C, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(C));
                                        removeCharsFromString(C, SI_chars);
                                        hfssBoundary.Capacitance = atof(
                                                                       C.data()) * factor;
                                    } else if (tmp_str.substr(0, 5)
                                               == "Faces") {
                                        hfssBoundary.faces.push_back(
                                            atoi(
                                                std::string(
                                                    tmp_str.substr(
                                                        6,
                                                        tmp_str.size()
                                                        - 7)).data()));
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            hfssBoundary.faces.push_back(
                                                atoi(
                                                    std::string(
                                                        tmp_str.substr(
                                                            0,
                                                            tmp_str.size()
                                                            - 1)).data()));
                                        }
                                    } else if (tmp_str.substr(0, 7)
                                               == "Objects") {
                                        hfssBoundary.solids.push_back(
                                            atoi(
                                                std::string(
                                                    tmp_str.substr(
                                                        8,
                                                        tmp_str.size()
                                                        - 9)).data()));
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            hfssBoundary.solids.push_back(
                                                atoi(
                                                    std::string(
                                                        tmp_str.substr(
                                                            0,
                                                            tmp_str.size()
                                                            - 1)).data()));
                                        }
                                    } else if (tmp_str == "$end") {
                                        iss >> tmp_str;
                                        if (tmp_str
                                                == "\'" + hfssBoundary.name
                                                + "\'") {
                                            bnds.push_back(hfssBoundary);
                                            break;
                                        }
                                    }
                                }
                            } else if (tmp_str == "$end") {
                                iss >> tmp_str;
                                if (tmp_str == "\'Boundaries\'") {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "Reading points";
// read_points();
    {
        size_t numPnts;
        int tmp_int;
        // double tmp_dbl, x, y, z;
        std::string tmp_str;
        std::string line;
        std::ifstream file_name(std::string(aux_path + "current.pnt").c_str());
        if (file_name.is_open()) {
            getline(file_name, line);
            std::istringstream iss(line);
            iss >> tmp_str;
            if (tmp_str == "points") {
                iss >> numPnts;
            } else {
                throw std::string("current.pnt is not of points type");
            }
            msh.n_nodes = numPnts;
            msh.nod_pos.resize(numPnts, std::vector<double>(3));
            for (size_t i = 0; i < numPnts; i++) {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                iss >> tmp_int;
                iss >> msh.nod_pos[i][0];
                iss >> msh.nod_pos[i][1];
                iss >> msh.nod_pos[i][2];
            }
            file_name.close();
        }
    }
    std::cout << ", faces";
// read_faces();
    {
        size_t numFaces;
        int tmp_int, /* n0, n1, n2,*/t0, t1;
        std::vector<size_t> bndTag;
        size_t bndTagNum;
        // double tmp_dbl;
        std::string tmp_str;
        std::string line;
        std::ifstream file_name(std::string(aux_path + "current.fac").c_str());
        if (file_name.is_open()) {
            getline(file_name, line);
            std::istringstream iss(line);
            iss >> tmp_str;
            if (tmp_str == "faces_v2") {
                iss >> numFaces;
            } else {
                throw std::string("current.fac is not of faces_v2 type");
            }
            msh.n_faces = numFaces;
            msh.fac_nodes.resize(numFaces, std::vector<size_t>(3));
            fac_hfss_tag.resize(numFaces);
            msh.fac_lab.resize(numFaces);
            msh.fac_adj_tet.resize(numFaces);
            for (size_t i = 0; i < numFaces; i++) {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str; // f
                iss >> tmp_int; // id
                iss >> msh.fac_nodes[i][0];
                iss >> msh.fac_nodes[i][1];
                iss >> msh.fac_nodes[i][2];
                iss >> tmp_str; // h
                iss >> t0;
                iss >> t1;
                iss >> bndTagNum;
                bndTag.resize(bndTagNum);
                for (size_t ibnd = 0; ibnd < bndTagNum; ibnd++) {
                    iss >> bndTag[ibnd];
                }
                fac_hfss_tag[i] = bndTag;
            }
            {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str != "end_face") {
                    throw std::string("end_face not found");
                }
            }
            {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "NumFaces") {
                    size_t tag, label;
                    iss >> bndTagNum;
                    for (size_t i = 0; i < bndTagNum; i++) {
                        getline(file_name, line);
                        std::istringstream iss(line);
                        iss >> tag;
                        iss >> label;
                        bnd_map[label].push_back(tag);
                    }
                }
            }
            file_name.close();
        }
    }
    std::cout << ", hydras";
// read_hydras();
    {
        size_t numHydras;
        int tmp_int, /* n0, n1, n2, n3, f0, f1, f2, f3,*/b0, b1, b2, b3;
        size_t l0, l1, l2, l3, l4, l5/*, s0*/;
        // double tmp_dbl;
        std::string tmp_str;
        std::string line;
        std::ifstream file_name(std::string(aux_path + "current.hyd").c_str());
        if (file_name.is_open()) {
            getline(file_name, line);
            std::istringstream iss(line);
            iss >> tmp_str;
            if (tmp_str == "hydras") {
                iss >> numHydras;
            } else {
                throw std::string("current.fac is not of faces_v2 type");
            }
            msh.n_tetras = numHydras;
            msh.tet_nodes.resize(numHydras, std::vector<size_t>(4));
            msh.tet_edges.resize(numHydras, std::vector<size_t>(6));
            msh.tet_faces.resize(numHydras, std::vector<size_t>(4));
            msh.tet_lab.resize(numHydras);
            hfss_id.resize(numHydras);
            for (size_t i = 0; i < numHydras; i++) {
                getline(file_name, line);
                {
                    std::istringstream iss(line);
                    iss >> tmp_str; // h
                    iss >> tmp_int; // id
                    iss >> msh.tet_nodes[i][0];
                    iss >> msh.tet_nodes[i][1];
                    iss >> msh.tet_nodes[i][2];
                    iss >> msh.tet_nodes[i][3];
                }
                getline(file_name, line);
                {
                    std::istringstream iss(line);
                    iss >> tmp_str; // f
                    iss >> msh.tet_faces[i][0];
                    iss >> msh.tet_faces[i][1];
                    iss >> msh.tet_faces[i][2];
                    iss >> msh.tet_faces[i][3];
                }
                getline(file_name, line);
                {
                    std::istringstream iss(line);
                    iss >> tmp_str; // b
                    iss >> tmp_int;
                    iss >> b0;
                    iss >> b1;
                    iss >> b2;
                    iss >> b3;
                }
                getline(file_name, line);
                {
                    std::istringstream iss(line);
                    iss >> tmp_str; // l
                    iss >> l0;
                    iss >> l1;
                    iss >> l2;
                    iss >> l3;
                    iss >> l4;
                    iss >> l5;
                }
                getline(file_name, line);
                {
                    std::istringstream iss(line);
                    iss >> tmp_str; // s
                    iss >> hfss_id[i];
                }
            }
            {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str != "end_hydra") {
                    throw std::string("end_hydra not found");
                }
            }
            file_name.close();
            for (size_t i = 0; i < msh.n_faces; i++) {
                for (unsigned int j = 0; j < 3; j++)
                    msh.fac_nodes[i][j] -= 1;
            }
            for (size_t i = 0; i < msh.n_tetras; i++) {
                for (unsigned int j = 0; j < 4; j++) {
                    msh.tet_nodes[i][j] -= 1;
                    msh.tet_faces[i][j] -= 1;
                }
            }
        }
    }
    std::cout << "\n";
// finalize_mesh();
    {
        tet_flag = std::vector<bool>(msh.n_tetras, false);
        fac_flag = std::vector<bool>(msh.n_faces, false);
        node_flag = std::vector<bool>(msh.n_nodes, false);
        size_t idx = 0, nidx = 0, /* eidx = 0,*/fidx = 0, tidx = 0;
        std::map<std::string, std::vector<size_t> > gSolidTets;
        for (std::vector<hfss_part>::iterator it = parts.begin();
                it != parts.end(); it++) {
            if (it->solve_inside) {
                mdl_mtrl mtr;
                mtr.name = it->name;
                mtr.type = it->material;
                mtr.epsr = mtrls[it->material].permittivity;
                mtr.mur = mtrls[it->material].permeability;
                mtr.sigma = mtrls[it->material].conductivity;
                mtr.tand = mtrls[it->material].dielectric_loss_tangent;
                for (size_t tid = 0; tid < msh.n_tetras; tid++) {
                    if (it->id == hfss_id[tid]) {
                        tet_flag[tid] = true;
                        mtr.label = frm.mtrls.size() + 1;
                        gSolidTets[it->name].push_back(tid);
                        msh.max_edg_marker = std::max(msh.max_edg_marker,
                                                      mtr.label);
                    }
                }
                if (gSolidTets[it->name].size() == 0) {
                    gSolidTets.erase(it->name);
                } else {
                    mtr.tetras.resize(gSolidTets[it->name].size());
                    size_t idx = 0;
                    for (std::vector<size_t>::iterator iter =
                                gSolidTets[it->name].begin();
                            iter != gSolidTets[it->name].end(); iter++) {
                        mtr.tetras[idx++] = *iter;
                        msh.tet_lab[*iter] = mtr.label;
                    }
                    std::cout << mtr.name << " " << mtr.type << "\n";
                    frm.mtrls.push_back(mtr);
                }
            }
        }
        for (std::vector<hfss_bnd>::iterator it = bnds.begin();
                it != bnds.end(); it++) {
            mdl_bc bc;
            bc.type = it->type;
            bc.name = it->name;
            bc.label = frm.bcs.size() + 1;
            std::cout << bc.name << " " << bc.type;
            if (strcmp(bc.type.data(), "WavePort") == 0) {
                bc.num_modes = it->num_modes;
                std::cout << " " << bc.num_modes;
            } else if (strcmp(bc.type.data(), "LumpedPort") == 0) {
                bc.lumped_impedance = std::complex<double>(it->FullResistance,
                                      it->FullReactance);
                std::cout << " " << bc.lumped_impedance;
            } else if (strcmp(bc.type.data(), "Impedance") == 0) {
                bc.surf_impedance = std::complex<double>(it->Resistance,
                                    it->Reactance);
                std::cout << " " << bc.surf_impedance;
            } else if (strcmp(bc.type.data(), "LumpedRLC") == 0) {
                bc.lumped_impedance = std::complex<double>(0.0, 0.0);
                bc.R = it->Resistance;
                bc.L = it->Inductance;
                bc.C = it->Capacitance;
                std::cout << " " << bc.R << " " << bc.L << " " << bc.C;
            }
            std::cout << "\n";
            frm.bcs.push_back(bc);
        }
        if (debug) {
            std::cout << "TetMap\n";
        }
        std::vector<size_t> tetMap(msh.n_tetras, UINT_MAX);
        tidx = 0;
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            if (tet_flag[tid]) {
                tetMap[tid] = tidx++;
            }
        }
        if (debug) {
            std::cout << "newTetNodes newTetFaces newTetLab " << tidx << "\n";
        }
        std::vector<std::vector<size_t> > newTetNodes(tidx,
                std::vector<size_t>(4));
        std::vector<std::vector<size_t> > newTetFaces(tidx,
                std::vector<size_t>(4));
        std::vector<int> newTetLab(tidx);
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            if (tet_flag[tid]) {
                newTetNodes[tetMap[tid]] = msh.tet_nodes[tid];
                newTetFaces[tetMap[tid]] = msh.tet_faces[tid];
                newTetLab[tetMap[tid]] = msh.tet_lab[tid];
            }
        }
        msh.n_tetras = tidx;
        msh.tet_nodes = newTetNodes;
        msh.tet_faces = newTetFaces;
        msh.tet_lab = newTetLab;
        if (debug) {
            std::cout << "tet_mtrl newTetras\n";
        }
        for (size_t mtrid = 0; mtrid < frm.mtrls.size(); mtrid++) {
            std::vector<size_t> newTetras;
            for (size_t tid = 0; tid < frm.mtrls[mtrid].tetras.size(); tid++) {
                if (tetMap[frm.mtrls[mtrid].tetras[tid]] < UINT_MAX) {
                    newTetras.push_back(tetMap[frm.mtrls[mtrid].tetras[tid]]);
                }
            }
            frm.mtrls[mtrid].tetras.clear();
            frm.mtrls[mtrid].tetras.resize(newTetras.size());
            for (size_t tid = 0; tid < newTetras.size(); tid++) {
                frm.mtrls[mtrid].tetras[tid] = newTetras[tid];
            }
        }
        // reorder nodes and faces
        if (debug) {
            std::cout << "nodMap " << msh.n_nodes << " facMap " << msh.n_faces
                      << "\n";
        }
        std::vector<size_t> nodMap(msh.n_nodes, UINT_MAX);
        std::vector<size_t> facMap(msh.n_faces, UINT_MAX);
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            for (size_t i = 0; i < 4; i++) {
                if (node_flag[msh.tet_nodes[tid][i]] == false) {
                    node_flag[msh.tet_nodes[tid][i]] = true;
                    nodMap[msh.tet_nodes[tid][i]] = nidx++;
                }
                if (fac_flag[msh.tet_faces[tid][i]] == false) {
                    fac_flag[msh.tet_faces[tid][i]] = true;
                    facMap[msh.tet_faces[tid][i]] = fidx++;
                }
            }
        }
        // finishing with nodes
        if (debug) {
            std::cout << "tet_nodes tet_faces\n";
        }
        for (size_t tid = 0; tid < msh.n_tetras; tid++) {
            for (size_t i = 0; i < 4; i++) {
                msh.tet_nodes[tid][i] = nodMap[msh.tet_nodes[tid][i]];
                msh.tet_faces[tid][i] = facMap[msh.tet_faces[tid][i]];
            }
            std::sort(msh.tet_nodes[tid].begin(), msh.tet_nodes[tid].end());
        }
        if (debug) {
            std::cout << "fac_nodes\n";
        }
        for (size_t fid = 0; fid < msh.n_faces; fid++) {
            for (size_t i = 0; i < 3; i++) {
                msh.fac_nodes[fid][i] = nodMap[msh.fac_nodes[fid][i]];
            }
            std::sort(msh.fac_nodes[fid].begin(), msh.fac_nodes[fid].end());
        }
        if (debug) {
            std::cout << "newNodPos " << nidx << "\n";
        }
        std::vector<std::vector<double> > newNodPos(nidx,
                std::vector<double>(3));
        for (size_t nid = 0; nid < msh.n_nodes; nid++) {
            if (node_flag[nid]) {
                newNodPos[nodMap[nid]] = msh.nod_pos[nid];
            }
        }
        if (debug) {
            std::cout << "newFacNodes " << fidx << "\n";
        }
        std::vector<std::vector<size_t> > newFacNodes(fidx,
                std::vector<size_t>(3));
        for (size_t fid = 0; fid < msh.n_faces; fid++) {
            if (fac_flag[fid]) {
                std::sort(msh.fac_nodes[fid].begin(), msh.fac_nodes[fid].end());
                newFacNodes[facMap[fid]] = msh.fac_nodes[fid];
            }
        }
        msh.nod_pos = newNodPos;
        msh.fac_nodes = newFacNodes;
        msh.n_nodes = nidx;
        msh.n_faces = fidx;
        // remain to assign face labels
        // reorder faces in tetrahedron
        // create edges
        if (debug) {
            std::cout << "fac_lab\n";
        }
        msh.fac_lab.assign(msh.n_faces, 0); // non boundary flag
        for (size_t fid = 0; fid < fac_hfss_tag.size(); fid++) {
            for (std::vector<hfss_bnd>::iterator it = bnds.begin();
                    it != bnds.end(); it++) {
                if (it->faces.size()) {
                    for (std::vector<size_t>::iterator itids =
                                it->faces.begin(); itids != it->faces.end();
                            itids++) {
                        if (std::find(fac_hfss_tag[fid].begin(),
                                      fac_hfss_tag[fid].end(), *itids)
                                != fac_hfss_tag[fid].end()) {
                            for (std::vector<mdl_bc>::iterator bndit =
                                        frm.bcs.begin(); bndit != frm.bcs.end();
                                    bndit++) {
                                if (bndit->name == it->name) {
                                    if (debug) {
                                        std::cout << " " << it->name << " "
                                                  << facMap[fid] << " ";
                                    }
                                    msh.fac_lab[facMap[fid]] = bndit->label;
                                    bndit->faces.push_back(facMap[fid]);
                                    msh.max_fac_marker = std::max(
                                                             msh.max_fac_marker, bndit->label);
                                }
                            }
                        }
                    }
                }
                if (debug) {
                    std::cout << " " << it->name << " " << facMap[fid] << " ";
                }
                // solid based boundaries
                if (it->solids.size()) {
                    for (std::vector<size_t>::iterator itids =
                                it->solids.begin(); itids != it->solids.end();
                            itids++) {
                        std::vector<size_t> cid = bnd_map[*itids];
                        for (size_t idd = 0; idd < cid.size(); idd++) {
                            if (std::find(fac_hfss_tag[fid].begin(),
                                          fac_hfss_tag[fid].end(), cid[idd])
                                    != fac_hfss_tag[fid].end()) {
                                for (std::vector<mdl_bc>::iterator bndit =
                                            frm.bcs.begin(); bndit != frm.bcs.end();
                                        bndit++) {
                                    if (bndit->name == it->name) {
                                        if (debug) {
                                            std::cout << " " << it->name << " "
                                                      << facMap[fid] << " ";
                                        }
                                        msh.fac_lab[facMap[fid]] = bndit->label;
                                        bndit->faces.push_back(facMap[fid]);
                                        msh.max_fac_marker = std::max(
                                                                 msh.max_fac_marker,
                                                                 bndit->label);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        /// reorder faces and create edges
        if (debug) {
            std::cout << "\ntet_faces fac_adj_tet\n";
        }
        std::vector<size_t> nodes(4), faces(4);
        fidx = 0;
        for (size_t tit = 0; tit < msh.n_tetras; tit++) {
            nodes[0] = msh.tet_nodes[tit][0];
            nodes[1] = msh.tet_nodes[tit][1];
            nodes[2] = msh.tet_nodes[tit][2];
            nodes[3] = msh.tet_nodes[tit][3];
            faces[0] = msh.tet_faces[tit][0];
            faces[1] = msh.tet_faces[tit][1];
            faces[2] = msh.tet_faces[tit][2];
            faces[3] = msh.tet_faces[tit][3];
            for (size_t ii = 0; ii < 4; ii++) {
                size_t fid = faces[ii];
                if ((nodes[1] == msh.fac_nodes[fid][0])
                        & (nodes[2] == msh.fac_nodes[fid][1])
                        & (nodes[3] == msh.fac_nodes[fid][2])) {
                    msh.tet_faces[tit][0] = fid;
                } else if ((nodes[0] == msh.fac_nodes[fid][0])
                           & (nodes[2] == msh.fac_nodes[fid][1])
                           & (nodes[3] == msh.fac_nodes[fid][2])) {
                    msh.tet_faces[tit][1] = fid;
                } else if ((nodes[0] == msh.fac_nodes[fid][0])
                           & (nodes[1] == msh.fac_nodes[fid][1])
                           & (nodes[3] == msh.fac_nodes[fid][2])) {
                    msh.tet_faces[tit][2] = fid;
                } else if ((nodes[0] == msh.fac_nodes[fid][0])
                           & (nodes[1] == msh.fac_nodes[fid][1])
                           & (nodes[2] == msh.fac_nodes[fid][2])) {
                    msh.tet_faces[tit][3] = fid;
                }
                msh.fac_adj_tet[fid].push_back(tit);
            }
        }
        /// create edges
        if (debug) {
            std::cout << "edgesMap\n";
        }
        idx = 0;
        size_t n0, n1;
        std::map<std::pair<size_t, size_t>, size_t> edgesMap;
        for (size_t tit = 0; tit < msh.n_tetras; tit++) {
            n0 = msh.fac_nodes[msh.tet_faces[tit][2]][1];
            n1 = msh.fac_nodes[msh.tet_faces[tit][2]][2];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = idx++;
            }
            n0 = msh.fac_nodes[msh.tet_faces[tit][1]][0];
            n1 = msh.fac_nodes[msh.tet_faces[tit][1]][2];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = idx++;
            }
            n0 = msh.fac_nodes[msh.tet_faces[tit][1]][1];
            n1 = msh.fac_nodes[msh.tet_faces[tit][1]][2];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = idx++;
            }
            n0 = msh.fac_nodes[msh.tet_faces[tit][0]][0];
            n1 = msh.fac_nodes[msh.tet_faces[tit][0]][1];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = idx++;
            }
            n0 = msh.fac_nodes[msh.tet_faces[tit][1]][0];
            n1 = msh.fac_nodes[msh.tet_faces[tit][1]][1];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = idx++;
            }
            n0 = msh.fac_nodes[msh.tet_faces[tit][2]][0];
            n1 = msh.fac_nodes[msh.tet_faces[tit][2]][1];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = idx++;
            }
        }
        msh.n_edges = edgesMap.size();
        if (debug) {
            std::cout << "edg_nodes\n";
        }
        std::map<std::pair<size_t, size_t>, size_t>::iterator emIter;
        msh.edg_nodes.resize(msh.n_edges, std::vector<size_t>(2));
        msh.edg_lab.resize(msh.n_edges, 0);
        for (emIter = edgesMap.begin(); emIter != edgesMap.end(); emIter++) {
            msh.edg_nodes[emIter->second][0] = emIter->first.first;
            msh.edg_nodes[emIter->second][1] = emIter->first.second;
        }
        idx = 0;
        if (debug) {
            std::cout << "fac_edges\n";
        }
        msh.fac_edges.resize(msh.n_faces, std::vector<size_t>(3));
        for (size_t fit = 0; fit < msh.n_faces; fit++) {
            n0 = msh.fac_nodes[fit][1];
            n1 = msh.fac_nodes[fit][2];
            msh.fac_edges[fit][0] = edgesMap[std::make_pair(n0, n1)];
            n0 = msh.fac_nodes[fit][0];
            n1 = msh.fac_nodes[fit][2];
            msh.fac_edges[fit][1] = edgesMap[std::make_pair(n0, n1)];
            n0 = msh.fac_nodes[fit][0];
            n1 = msh.fac_nodes[fit][1];
            msh.fac_edges[fit][2] = edgesMap[std::make_pair(n0, n1)];
        }
        if (debug) {
            std::cout << "fac_edges tet_edges\n";
        }
        for (size_t tit = 0; tit < msh.n_tetras; tit++) {
            msh.tet_edges[tit][0] = msh.fac_edges[msh.tet_faces[tit][2]][2];
            msh.tet_edges[tit][1] = msh.fac_edges[msh.tet_faces[tit][1]][2];
            msh.tet_edges[tit][2] = msh.fac_edges[msh.tet_faces[tit][1]][1];
            msh.tet_edges[tit][3] = msh.fac_edges[msh.tet_faces[tit][0]][2];
            msh.tet_edges[tit][4] = msh.fac_edges[msh.tet_faces[tit][0]][1];
            msh.tet_edges[tit][5] = msh.fac_edges[msh.tet_faces[tit][0]][0];
        }
    }
}

void mdl_core::wrap_aedt(std::string& name, std::string& aux_path) {
    struct aedt_part {
        std::string name;
        std::string material;
        bool solve_inside = false;
        size_t id;
    };
    struct aedt_bnd {
        std::string name;
        std::string type;
        std::vector<size_t> faces;
        std::vector<size_t> solids;
        std::vector<size_t> face_ids;
        int num_modes;
        double FullResistance = 0.0;
        double FullReactance = 0.0;
        double Resistance = 0.0;
        double Reactance = 0.0;
        double Inductance = 0.0;
        double Capacitance = 0.0;
    };
    struct aedt_mtrl {
        double permittivity = 1.0;
        double permeability = 1.0;
        double conductivity = 0.0;
        double dielectric_loss_tangent = 0.0;
        std::string name;
    };
    bool debug = false;
    std::map<std::string, aedt_mtrl> mtrls;
    std::vector<size_t> mtrl_tag;
    std::vector<size_t> aedt_id;
    std::vector<bool> tet_flag;
    std::vector<bool> fac_flag;
    std::vector<bool> node_flag;
    std::vector<std::vector<size_t> > fac_aedt_tag;
    std::vector<aedt_bnd> bnds;
    std::vector<aedt_part> parts;
    std::map<size_t, std::vector<size_t> > bnd_map;
    std::vector<std::vector<size_t> > adj_tetra;

// starting algorithm
    msh.clear();
    msh.type = "TETRA";
    frm.type = "EM_E_FD";
    std::cout << "- Parsing AEDT project files\n";
//  READ MAIN aedt
    {
        std::string R, L, C, X;
        double tmp_dbl;
        std::string tmp_str;
        std::string partName;
        // bool solve_inside;
        std::string materialName;
        std::string boundaryName;
        std::string line;
        std::ifstream file_name(std::string(name + ".aedt").c_str());
        if (file_name.is_open()) {
            while (getline(file_name, line)) { //file_name.good())
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "$begin") {
                    iss >> tmp_str;
                    if (tmp_str == "\'Materials\'") {
                        while (getline(file_name, line)) {
                            std::istringstream iss(line);
                            iss >> tmp_str;
                            if (tmp_str == "$begin") {
                                aedt_mtrl aedtMaterial;
                                std::string tmpMtrl;
                                while (iss.good()) {
                                    iss >> tmp_str;
                                    tmpMtrl.append(tmp_str);
                                }
                                materialName = tmpMtrl.substr(1,
                                                              tmpMtrl.size() - 2);
                                while (getline(file_name, line)) {
                                    std::istringstream iss(line);
                                    iss >> tmp_str; // $begin
                                    if (tmp_str.substr(0, 12)
                                            == "permittivity") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        14,
                                                        tmp_str.size()
                                                        - 15)).data());
                                        aedtMaterial.permittivity = tmp_dbl;
                                    } else if (tmp_str.substr(0, 12)
                                               == "permeability") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        14,
                                                        tmp_str.size()
                                                        - 15)).data());
                                        aedtMaterial.permeability = tmp_dbl;
                                    } else if (tmp_str.substr(0, 12)
                                               == "conductivity") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        14,
                                                        tmp_str.size()
                                                        - 15)).data());
                                        aedtMaterial.conductivity = tmp_dbl;
                                    } else if (tmp_str.substr(0, 23)
                                               == "dielectric_loss_tangent") {
                                        tmp_dbl =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        25,
                                                        tmp_str.size()
                                                        - 26)).data());
                                        aedtMaterial.dielectric_loss_tangent =
                                            tmp_dbl;
                                    } else if (tmp_str == "$end") {
                                        std::string tmpMtrl;
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            tmpMtrl.append(tmp_str);
                                        }
                                        tmpMtrl = tmpMtrl.substr(1,
                                                                 tmpMtrl.size() - 2);
                                        if (tmpMtrl == materialName) {
                                            aedtMaterial.name = materialName;
                                            mtrls[materialName] = aedtMaterial;
                                            break;
                                        }
                                    }
                                }
                            } else if (tmp_str == "$end") {
                                iss >> tmp_str;
                                if (tmp_str == "\'Materials\'") {
                                    break;
                                }
                            }
                        }
                    } else if (tmp_str == "\'ToplevelParts\'") {
                        while (getline(file_name, line)) {
                            std::istringstream iss(line);
                            iss >> tmp_str;
                            if (tmp_str == "$begin") {
                                iss >> tmp_str;
                                if (tmp_str == "\'GeometryPart\'") {
                                    aedt_part aedtPart;
                                    while (getline(file_name, line)) {
                                        std::istringstream iss(line);
                                        iss >> tmp_str;
                                        if (tmp_str == "$begin") {
                                            iss >> tmp_str;
                                            if (tmp_str == "\'Attributes\'") {
                                                while (getline(file_name, line)) {
                                                    std::istringstream iss(
                                                        line);
                                                    iss >> tmp_str;
                                                    if (tmp_str.substr(0, 4)
                                                            == "Name") {
                                                        std::string tmpString =
                                                            tmp_str;
                                                        while (iss.good()) {
                                                            iss >> tmp_str;
                                                            tmpString.append(
                                                                tmp_str);
                                                        }
                                                        aedtPart.name =
                                                            std::string(
                                                                tmpString.substr(
                                                                    6,
                                                                    tmpString.size()
                                                                    - 7));
                                                    } else if (tmp_str.substr(0,
                                                                              13)
                                                               == "MaterialValue") { // for hfss v13
                                                        std::string tmpString =
                                                            tmp_str;
                                                        while (iss.good()) {
                                                            iss >> tmp_str;
                                                            tmpString.append(
                                                                tmp_str);
                                                        }
                                                        aedtPart.material =
                                                            std::string(
                                                                tmpString.substr(
                                                                    16,
                                                                    tmpString.size()
                                                                    - 18).data());
                                                    } else if (tmp_str.substr(0,
                                                                              12)
                                                               == "MaterialName") { // for hfss v11
                                                        std::string tmpString =
                                                            tmp_str;
                                                        while (iss.good()) {
                                                            iss >> tmp_str;
                                                            tmpString.append(
                                                                tmp_str);
                                                        }
                                                        aedtPart.material =
                                                            std::string(
                                                                tmpString.substr(
                                                                    14,
                                                                    tmpString.size()
                                                                    - 15));
                                                    } else if (tmp_str.substr(0,
                                                                              11)
                                                               == "SolveInside") {
                                                        aedtPart.solve_inside =
                                                            std::string(
                                                                tmp_str.substr(
                                                                    12,
                                                                    tmp_str.size()
                                                                    - 12))
                                                            == "true";
                                                    } else if (tmp_str
                                                               == "$end") {
                                                        iss >> tmp_str;
                                                        if (tmp_str
                                                                == "\'Attributes\'") {
                                                            break;
                                                        }
                                                    }
                                                }
                                            } else if (tmp_str
                                                       == "\'Operation\'") {
                                                while (getline(file_name, line)) {
                                                    std::istringstream iss(
                                                        line);
                                                    iss >> tmp_str;
                                                    if (tmp_str.substr(0, 12)
                                                            == "ParentPartID") {
                                                        aedtPart.id =
                                                            atoi(
                                                                std::string(
                                                                    tmp_str.substr(
                                                                        13,
                                                                        tmp_str.size()
                                                                        - 13)).data());
                                                    } else if (tmp_str
                                                               == "$end") {
                                                        iss >> tmp_str;
                                                        if (tmp_str
                                                                == "\'Operation\'") {
                                                            break;
                                                        }
                                                    }
                                                }
                                            }
                                        } else if (tmp_str == "$end") {
                                            iss >> tmp_str;
                                            if (tmp_str == "\'GeometryPart\'") {
                                                parts.push_back(aedtPart);
                                                break;
                                            }
                                        }
                                    }
                                }
                            } else if (tmp_str == "$end") {
                                iss >> tmp_str;
                                if (tmp_str == "\'ToplevelParts\'") {
                                    break;
                                }
                            }
                        }
                    } else if (tmp_str == "\'Boundaries\'") {
                        while (getline(file_name, line)) {
                            std::istringstream iss(line);
                            iss >> tmp_str; // $begin
                            if (tmp_str == "$begin") {
                                iss >> tmp_str;
                                aedt_bnd aedtBoundary;
                                aedtBoundary.name = tmp_str.substr(1,
                                                                   tmp_str.size() - 2);
                                while (getline(file_name, line)) {
                                    std::istringstream iss(line);
                                    iss >> tmp_str; // $begin
                                    if (tmp_str.substr(0, 9) == "BoundType") {
                                        aedtBoundary.type = std::string(
                                                                tmp_str.substr(11,
                                                                               tmp_str.size() - 11));
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            aedtBoundary.type += std::string(
                                                                     tmp_str.substr(0,
                                                                                    tmp_str.size()
                                                                                    - 1));
                                        }
                                        if (aedtBoundary.type.substr(0, 9)
                                                == "Radiation") {
                                            aedtBoundary.type = "Radiation";
                                        } else if (aedtBoundary.type.substr(0,
                                                                            9) == "Impedance") {
                                            aedtBoundary.type = "Impedance";
                                        } else if (aedtBoundary.type.substr(0,
                                                                            7) == "Voltage") {
                                            aedtBoundary.type = "Voltage";
                                        } else if (aedtBoundary.type.substr(0,
                                                                            8) == "Terminal") {
                                            aedtBoundary.type = "Terminal";
                                        }
                                    } else if (tmp_str.substr(0, 8)
                                               == "NumModes") {
                                        aedtBoundary.num_modes =
                                            atoi(
                                                std::string(
                                                    tmp_str.substr(
                                                        9,
                                                        tmp_str.size()
                                                        - 1)).data());
                                    } else if (tmp_str.substr(0, 14)
                                               == "FullResistance") {
                                        aedtBoundary.FullResistance =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        16,
                                                        tmp_str.size()
                                                        - 20)).data());
                                    } else if (tmp_str.substr(0, 13)
                                               == "FullReactance") {
                                        aedtBoundary.FullReactance =
                                            atof(
                                                std::string(
                                                    tmp_str.substr(
                                                        15,
                                                        tmp_str.size()
                                                        - 19)).data());
                                    } else if (tmp_str.substr(0, 10)
                                               == "Resistance") {
                                        R = tmp_str.substr(11,
                                                           tmp_str.size() - 11).data();
                                        removeCharsFromString(R, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(R));
                                        removeCharsFromString(R, SI_chars);
                                        aedtBoundary.Resistance = atof(R.data())
                                                                  * factor;
                                    } else if (tmp_str.substr(0, 9)
                                               == "Reactance") {
                                        X = tmp_str.substr(10,
                                                           tmp_str.size() - 10).data();
                                        removeCharsFromString(X, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(X));
                                        removeCharsFromString(X, SI_chars);
                                        aedtBoundary.Reactance = atof(X.data())
                                                                 * factor;
                                    } else if (tmp_str.substr(0, 10)
                                               == "Inductance") {
                                        L = tmp_str.substr(11,
                                                           tmp_str.size() - 11).data();
                                        removeCharsFromString(L, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(L));
                                        removeCharsFromString(L, SI_chars);
                                        aedtBoundary.Inductance = atof(L.data())
                                                                  * factor;
                                    } else if (tmp_str.substr(0, 11)
                                               == "Capacitance") {
                                        C = tmp_str.substr(12,
                                                           tmp_str.size() - 12).data();
                                        removeCharsFromString(C, "\'");
                                        double factor = set_factor(
                                                            find_SI_factor(C));
                                        removeCharsFromString(C, SI_chars);
                                        aedtBoundary.Capacitance = atof(
                                                                       C.data()) * factor;
                                    } else if (tmp_str.substr(0, 5)
                                               == "Faces") {
                                        aedtBoundary.faces.push_back(
                                            atoi(
                                                std::string(
                                                    tmp_str.substr(
                                                        6,
                                                        tmp_str.size()
                                                        - 7)).data()));
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            aedtBoundary.faces.push_back(
                                                atoi(
                                                    std::string(
                                                        tmp_str.substr(
                                                            0,
                                                            tmp_str.size()
                                                            - 1)).data()));
                                        }
                                    } else if (tmp_str.substr(0, 7)
                                               == "Objects") {
                                        aedtBoundary.solids.push_back(
                                            atoi(
                                                std::string(
                                                    tmp_str.substr(
                                                        8,
                                                        tmp_str.size()
                                                        - 9)).data()));
                                        while (iss.good()) {
                                            iss >> tmp_str;
                                            aedtBoundary.solids.push_back(
                                                atoi(
                                                    std::string(
                                                        tmp_str.substr(
                                                            0,
                                                            tmp_str.size()
                                                            - 1)).data()));
                                        }
                                    } else if (tmp_str == "$end") {
                                        iss >> tmp_str;
                                        if (tmp_str
                                                == "\'" + aedtBoundary.name
                                                + "\'") {
                                            bnds.push_back(aedtBoundary);
                                            break;
                                        }
                                    }
                                }
                            } else if (tmp_str == "$end") {
                                iss >> tmp_str;
                                if (tmp_str == "\'Boundaries\'") {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    {
        // read ngmesh and process model data
        std::map<int, std::vector<int> > map_body_id_face_ids;
        std::map<int, int> map_facet_id_body_id, map_facet_id_face_id;
        int tmp_int;
        std::string tmp_str;
        std::string line;
        std::ifstream file_name(
            std::string(aux_path + "current.ngmesh").c_str());
        if (file_name.is_open()) {
            int body_id;
            do {
                getline(file_name, line);
            } while (line.compare(0, 15, "begin_body_data") != 0);
            std::cout << "Reading bodies";
            do {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "body_id") {
                    iss >> body_id;
                    getline(file_name, line);
                    iss.clear();
                    iss.str("");
                    iss.str(line);
                    iss >> tmp_str;
                    //std::cout << "\n" << body_id << ":";
                    if (tmp_str == "face_ids")
                        while (iss >> tmp_int) {
                            map_body_id_face_ids[body_id].push_back(tmp_int);
                            //std::cout << " " << tmp_int;
                        }
                }
            } while (line.compare(0, 13, "end_body_data") != 0);
            do {
                getline(file_name, line);
            } while (line.compare(0, 16, "begin_point_data") != 0);
            std::cout << ", points";
            do {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "pid") {
                    iss >> tmp_int; // id
                    iss >> tmp_str; // coords
                    std::vector<double> tmp_coord(3);
                    iss >> tmp_coord[0];
                    iss >> tmp_coord[1];
                    iss >> tmp_coord[2];
                    msh.nod_pos.push_back(tmp_coord);
                    msh.n_nodes++;
                }
            } while (line.compare(0, 14, "end_point_data") != 0);
            if (debug)
                std::cout << "\nn_nodes = " << msh.n_nodes << "\n";
            do {
                getline(file_name, line);
            } while (line.compare(0, 16, "begin_facet_data") != 0);
            std::cout << ", facets";

            int facet_id;
            do {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "facet_id") {
                    iss >> facet_id;
                    for (int i = 0; i < 13; i++) {
                        iss >> tmp_str;
                    }
                    iss >> tmp_int;
                    map_facet_id_face_id[facet_id] = tmp_int;
                    for (std::map<int, std::vector<int> >::iterator it =
                                map_body_id_face_ids.begin();
                            it != map_body_id_face_ids.end(); it++)
                        if (std::find(it->second.begin(), it->second.end(),
                                      tmp_int) != it->second.end()) {
                            map_facet_id_body_id[facet_id] = it->first;
                        }
                }
                if (tmp_str == "seid") {
                    iss >> tmp_int; // id
                    iss >> tmp_str; // facet_id
                    int tmp_facet_id;
                    iss >> tmp_facet_id;
                    msh.fac_lab.push_back(map_facet_id_face_id[tmp_facet_id]);
                    iss >> tmp_str; // nv
                    iss >> tmp_int; // 3
                    iss >> tmp_str; // vert_ids
                    std::vector<size_t> tmp_facet(3);
                    iss >> tmp_facet[0];
                    iss >> tmp_facet[1];
                    iss >> tmp_facet[2];
                    sort(tmp_facet.begin(), tmp_facet.end());
                    msh.fac_nodes.push_back(tmp_facet);
                    msh.n_faces++;
                }
            } while (line.compare(0, 14, "end_facet_data") != 0);
            if (debug)
                std::cout << "\nn_faces = " << msh.n_faces << "\n";
            msh.fac_lab.resize(msh.n_faces);
            msh.fac_adj_tet.resize(msh.n_faces);
            do {
                getline(file_name, line);
            } while (line.compare(0, 22, "begin_vol_element_data") != 0);
            std::cout << ", vols";
            do {
                getline(file_name, line);
                std::istringstream iss(line);
                iss >> tmp_str;
                if (tmp_str == "veid") {
                    iss >> tmp_int; // id
                    iss >> tmp_str; // body_id
                    int tmp_body_id;
                    iss >> tmp_body_id;
                    msh.tet_lab.push_back(tmp_body_id);
                    iss >> tmp_str; // nv
                    iss >> tmp_int; // 4
                    iss >> tmp_str; // vert_ids
                    std::vector<size_t> tmp_vol(4);
                    iss >> tmp_vol[0];
                    iss >> tmp_vol[1];
                    iss >> tmp_vol[2];
                    iss >> tmp_vol[3];
                    sort(tmp_vol.begin(), tmp_vol.end());
                    msh.tet_nodes.push_back(tmp_vol);
                    msh.n_tetras++;
                }
            } while (line.compare(0, 20, "end_vol_element_data") != 0);
            if (debug)
                std::cout << "\nn_tetras = " << msh.n_tetras << "\n";
            msh.tet_edges.resize(msh.n_tetras, std::vector<size_t>(6));
            msh.tet_faces.resize(msh.n_tetras, std::vector<size_t>(4));
            msh.tet_lab.resize(msh.n_tetras);
            file_name.close();
        }
        for (size_t i = 0; i < msh.n_faces; i++) {
            for (unsigned int j = 0; j < 3; j++)
                msh.fac_nodes[i][j] -= 1;
        }
        for (size_t i = 0; i < msh.n_tetras; i++) {
            for (unsigned int j = 0; j < 4; j++) {
                msh.tet_nodes[i][j] -= 1;
            }
        }
        std::vector<int> bc_ids = msh.fac_lab, mtrl_ids = msh.tet_lab;
        bc_ids.erase(unique(bc_ids.begin(), bc_ids.end()), bc_ids.end());
        mtrl_ids.erase(unique(mtrl_ids.begin(), mtrl_ids.end()),
                       mtrl_ids.end());
        // for(int n : bc_ids) std::cout << n << " ";
        std::cout << "\n";
        tet_flag.assign(msh.n_tetras, false);
        for (std::vector<aedt_part>::iterator it = parts.begin();
                it != parts.end(); it++) {
            if (it->solve_inside) {
                mdl_mtrl mtr;
                mtr.name = it->name;
                mtr.type = it->material;
                mtr.epsr = mtrls[it->material].permittivity;
                mtr.mur = mtrls[it->material].permeability;
                mtr.sigma = mtrls[it->material].conductivity;
                mtr.tand = mtrls[it->material].dielectric_loss_tangent;
                for (size_t tid = 0; tid < msh.n_tetras; tid++) {
                    if (it->id == msh.tet_lab[tid]) {
                        tet_flag[tid] = true;
                        mtr.label = frm.mtrls.size() + 1;
                        msh.tet_lab[tid] = mtr.label;
                    }
                }
                if (std::find(mtrl_ids.begin(), mtrl_ids.end(), it->id)
                        != mtrl_ids.end()) {
                    std::cout << mtr.name << " " << mtr.type << " " << mtr.label
                              << "\n";
                    frm.mtrls.push_back(mtr);
                }
            } else {
                aedt_bnd bc;
                bc.name = it->name;
                bc.type = "PerfectE";
                bc.solids.push_back(it->id);
                bnds.push_back(bc);
            }
        }
        {
            // remove superfluous tets
            size_t cnt = 0;
            std::vector<std::vector<size_t> > new_tet_nodes;
            std::vector<int> new_tet_lab;
            for (size_t tid = 0; tid < msh.n_tetras; tid++) {
                if (tet_flag[tid]) {
                    new_tet_nodes.push_back(msh.tet_nodes[tid]);
                    new_tet_lab.push_back(msh.tet_lab[tid]);
                    cnt++;
                }
            }
            msh.n_tetras = cnt;
            msh.tet_nodes.swap(new_tet_nodes);
            msh.tet_lab.swap(new_tet_lab);
        }
        fac_flag.assign(msh.n_faces, true);
        for (std::vector<aedt_bnd>::iterator it = bnds.begin();
                it != bnds.end(); it++) {
            mdl_bc bc;
            bc.type = it->type;
            bc.name = it->name;
            bc.label = frm.bcs.size() + 1;
            std::cout << bc.name << " " << bc.type << " " << bc.label;
            for (int n : it->faces)
                for (size_t i = 0; i < msh.n_faces; i++)
                    if (fac_flag[i])
                        if (msh.fac_lab[i] == n) {
                            fac_flag[i] = false;
                            msh.fac_lab[i] = bc.label;
                        }
            for (int n : it->solids)
                for (size_t i = 0; i < msh.n_faces; i++)
                    if (fac_flag[i])
                        for (int j = 0; j < map_body_id_face_ids[n].size(); j++)
                            if (msh.fac_lab[i] == map_body_id_face_ids[n][j]) {
                                fac_flag[i] = false;
                                msh.fac_lab[i] = bc.label;
                            }
            if (strcmp(bc.type.data(), "WavePort") == 0) {
                bc.num_modes = it->num_modes;
                std::cout << " " << bc.num_modes;
            } else if (strcmp(bc.type.data(), "LumpedPort") == 0) {
                bc.lumped_impedance = std::complex<double>(it->FullResistance,
                                      it->FullReactance);
                std::cout << " " << bc.lumped_impedance;
            } else if (strcmp(bc.type.data(), "Impedance") == 0) {
                bc.surf_impedance = std::complex<double>(it->Resistance,
                                    it->Reactance);
                std::cout << " " << bc.surf_impedance;
            } else if (strcmp(bc.type.data(), "LumpedRLC") == 0) {
                bc.lumped_impedance = std::complex<double>(0.0, 0.0);
                bc.R = it->Resistance;
                bc.L = it->Inductance;
                bc.C = it->Capacitance;
                std::cout << " " << bc.R << " " << bc.L << " " << bc.C;
            }
            std::cout << "\n";
            frm.bcs.push_back(bc);
        }
        int id = frm.bcs.size() + 1;
        std::vector<int> face_ids;
        std::map<int, int> map_body_bc;
        for (size_t i = 0; i < msh.n_faces; i++)
            if (fac_flag[i]) {
                face_ids.push_back(msh.fac_lab[i]);
                face_ids.erase(unique(face_ids.begin(), face_ids.end()),
                               face_ids.end());
            }
        for (int n : face_ids) {
            for (std::map<int, std::vector<int> >::iterator it =
                        map_body_id_face_ids.begin();
                    it != map_body_id_face_ids.end(); it++)
                if (std::find(it->second.begin(), it->second.end(), n)
                        != it->second.end())
                    map_body_bc[it->first] = id++;
        }
        std::vector<int> new_bcs;
        bool add_bcs = false;
        for (size_t i = 0; i < msh.n_faces; i++)
            if (fac_flag[i]) {
                msh.fac_lab[i] =
                    map_body_bc[map_facet_id_body_id[msh.fac_lab[i]]];
                new_bcs.push_back(msh.fac_lab[i]);
                std::sort(new_bcs.begin(), new_bcs.end());
                new_bcs.erase(unique(new_bcs.begin(), new_bcs.end()),
                              new_bcs.end());
                add_bcs = true;
            }
        id = frm.bcs.size() + 1;
        for (int j = 0; j < new_bcs.size(); j++) {
            for (size_t i = 0; i < msh.n_faces; i++)
                if (fac_flag[i] && msh.fac_lab[i] == new_bcs[j])
                    msh.fac_lab[i] = id;
            new_bcs[j] = id++;
        }
        if (add_bcs) {
            for (int n : new_bcs) {
                mdl_bc bc;
                bc.type = "None";
                bc.name = "Bodies";
                bc.label = n;
                std::cout << bc.name << " " << bc.type << " " << bc.label
                          << "\n";
                frm.bcs.push_back(bc);
            }
        }
        {
            //removing superfluous nodes and faces
            std::vector<std::vector<double> > new_nod_pos;
            std::vector<bool> nod_flag(msh.n_nodes, false);
            std::vector<size_t> nod_map(msh.n_nodes);
            size_t rnid = 0;
            for (size_t tid = 0; tid < msh.n_tetras; tid++) {
                for (unsigned int i = 0; i < 4; i++) {
                    if (!nod_flag[msh.tet_nodes[tid][i]]) {
                        new_nod_pos.push_back(
                            msh.nod_pos[msh.tet_nodes[tid][i]]);
                        nod_map[msh.tet_nodes[tid][i]] = rnid++;
                        nod_flag[msh.tet_nodes[tid][i]] = true;
                    }
                }
            }
            msh.n_nodes = new_nod_pos.size();
            msh.nod_pos.swap(new_nod_pos);

            std::vector<std::vector<size_t> > new_fac_nodes;
            std::vector<int> new_fac_lab;
            for (size_t fid = 0; fid < msh.n_faces; fid++) {
                bool append = true;
                for (unsigned int i = 0; i < 3; i++)
                    if (!nod_flag[msh.fac_nodes[fid][i]])
                        append = false;
                if (append) {
                    new_fac_nodes.push_back(msh.fac_nodes[fid]);
                    new_fac_lab.push_back(msh.fac_lab[fid]);
                }
            }
            msh.n_faces = new_fac_nodes.size();
            msh.fac_nodes.swap(new_fac_nodes);
            msh.fac_lab.swap(new_fac_lab);

            for (size_t fid = 0; fid < msh.n_faces; fid++)
                for (unsigned int i = 0; i < 3; i++)
                    msh.fac_nodes[fid][i] = nod_map[msh.fac_nodes[fid][i]];
            for (size_t tid = 0; tid < msh.n_tetras; tid++)
                for (unsigned int i = 0; i < 4; i++)
                    msh.tet_nodes[tid][i] = nod_map[msh.tet_nodes[tid][i]];

        }

        // creating edges and internal faces
        size_t eidx = 0, fidx = 0;
        size_t n0, n1, n2;
        std::vector<size_t> tmp_edg(2), tmp_fac(3);
        std::map<std::pair<size_t, size_t>, size_t> edgesMap;
        std::map<std::tuple<size_t, size_t, size_t>, size_t> facesMap;
        for (size_t fid = 0; fid < msh.n_faces; fid++) {
            n0 = msh.fac_nodes[fid][0];
            n1 = msh.fac_nodes[fid][1];
            n2 = msh.fac_nodes[fid][2];
            facesMap[std::make_tuple(n0, n1, n2)] = fid;
        }
        fidx = msh.n_faces;
        for (size_t tit = 0; tit < msh.n_tetras; tit++) {
            // edges
            n0 = msh.tet_nodes[tit][0];
            n1 = msh.tet_nodes[tit][1];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = eidx++;
                tmp_edg[0] = n0;
                tmp_edg[1] = n1;
                msh.edg_nodes.push_back(tmp_edg);
            }
            n0 = msh.tet_nodes[tit][0];
            n1 = msh.tet_nodes[tit][2];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = eidx++;
                tmp_edg[0] = n0;
                tmp_edg[1] = n1;
                msh.edg_nodes.push_back(tmp_edg);
            }
            n0 = msh.tet_nodes[tit][0];
            n1 = msh.tet_nodes[tit][3];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = eidx++;
                tmp_edg[0] = n0;
                tmp_edg[1] = n1;
                msh.edg_nodes.push_back(tmp_edg);
            }
            n0 = msh.tet_nodes[tit][1];
            n1 = msh.tet_nodes[tit][2];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = eidx++;
                tmp_edg[0] = n0;
                tmp_edg[1] = n1;
                msh.edg_nodes.push_back(tmp_edg);
            }
            n0 = msh.tet_nodes[tit][1];
            n1 = msh.tet_nodes[tit][3];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = eidx++;
                tmp_edg[0] = n0;
                tmp_edg[1] = n1;
                msh.edg_nodes.push_back(tmp_edg);
            }
            n0 = msh.tet_nodes[tit][2];
            n1 = msh.tet_nodes[tit][3];
            if (edgesMap.find(std::make_pair(n0, n1)) == edgesMap.end()) {
                edgesMap[std::make_pair(n0, n1)] = eidx++;
                tmp_edg[0] = n0;
                tmp_edg[1] = n1;
                msh.edg_nodes.push_back(tmp_edg);
            }
            // faces
            n0 = msh.tet_nodes[tit][1];
            n1 = msh.tet_nodes[tit][2];
            n2 = msh.tet_nodes[tit][3];
            if (facesMap.find(std::make_tuple(n0, n1, n2)) == facesMap.end()) {
                facesMap[std::make_tuple(n0, n1, n2)] = fidx++;
                tmp_fac[0] = n0;
                tmp_fac[1] = n1;
                tmp_fac[2] = n2;
                msh.fac_nodes.push_back(tmp_fac);
                //msh.fac_lab.push_back(0);
            }
            n0 = msh.tet_nodes[tit][0];
            n1 = msh.tet_nodes[tit][2];
            n2 = msh.tet_nodes[tit][3];
            if (facesMap.find(std::make_tuple(n0, n1, n2)) == facesMap.end()) {
                facesMap[std::make_tuple(n0, n1, n2)] = fidx++;
                tmp_fac[0] = n0;
                tmp_fac[1] = n1;
                tmp_fac[2] = n2;
                msh.fac_nodes.push_back(tmp_fac);
                //msh.fac_lab.push_back(0);
            }
            n0 = msh.tet_nodes[tit][0];
            n1 = msh.tet_nodes[tit][1];
            n2 = msh.tet_nodes[tit][3];
            if (facesMap.find(std::make_tuple(n0, n1, n2)) == facesMap.end()) {
                facesMap[std::make_tuple(n0, n1, n2)] = fidx++;
                tmp_fac[0] = n0;
                tmp_fac[1] = n1;
                tmp_fac[2] = n2;
                msh.fac_nodes.push_back(tmp_fac);
                //msh.fac_lab.push_back(0);
            }
            n0 = msh.tet_nodes[tit][0];
            n1 = msh.tet_nodes[tit][1];
            n2 = msh.tet_nodes[tit][2];
            if (facesMap.find(std::make_tuple(n0, n1, n2)) == facesMap.end()) {
                facesMap[std::make_tuple(n0, n1, n2)] = fidx++;
                tmp_fac[0] = n0;
                tmp_fac[1] = n1;
                tmp_fac[2] = n2;
                msh.fac_nodes.push_back(tmp_fac);
                //msh.fac_lab.push_back(0);
            }
        }
        msh.n_faces = fidx;
        msh.n_edges = eidx;
        msh.fac_lab.resize(fidx);
        msh.edg_lab.resize(eidx);
        msh.regularize_mesh();
        // msh.get_mesh_statistics();
    }
}


void mdl_core::create_tri_mesh() {
    msh.clear();
    msh.type = "TRIA";
    msh.n_nodes = sld.nodes.size();
    msh.nod_pos = sld.nodes;
    msh.n_faces = sld.faces.size();
    msh.fac_nodes.resize(msh.n_faces);
    msh.fac_edges.resize(msh.n_faces);
    msh.fac_lab.assign(msh.n_faces,1);
    std::map<std::pair<size_t,size_t>, size_t> edg_map;
    size_t edg_cnt = 0;
    for(size_t i = 0; i < msh.n_faces; i++) {
        msh.fac_nodes[i] = sld.faces[i].polygons[0];
        std::sort(msh.fac_nodes[i].begin(), msh.fac_nodes[i].end());
//		msh.fac_lab[i] = sld.faces_marker[i][0];
        if(edg_map.find(std::make_pair(msh.fac_nodes[i][0],msh.fac_nodes[i][1]))
                == edg_map.end())
            edg_map[std::make_pair(msh.fac_nodes[i][0],msh.fac_nodes[i][1])] = edg_cnt++;
        if(edg_map.find(std::make_pair(msh.fac_nodes[i][0],msh.fac_nodes[i][1]))
                == edg_map.end())
            edg_map[std::make_pair(msh.fac_nodes[i][0],msh.fac_nodes[i][2])] = edg_cnt++;
        if(edg_map.find(std::make_pair(msh.fac_nodes[i][0],msh.fac_nodes[i][1]))
                == edg_map.end())
            edg_map[std::make_pair(msh.fac_nodes[i][1],msh.fac_nodes[i][2])] = edg_cnt++;
    }
    std::cout << edg_cnt << "\n";
    msh.n_edges = edg_cnt;
    msh.edg_nodes.resize(edg_cnt);
    for(std::map<std::pair<size_t,size_t>, size_t>::iterator it = edg_map.begin();
            it != edg_map.end(); it++) {
        std::vector<size_t> edge(2);
        edge[0] = std::get <0>(it->first);
        edge[1] = std::get <1>(it->first);
        msh.edg_nodes[it->second] = edge;
    }
    msh.edg_lab.assign(edg_cnt,1);
    msh.max_edg_marker = 1;
    msh.get_mesh_statistics();
}

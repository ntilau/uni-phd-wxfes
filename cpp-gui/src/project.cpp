#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>


#ifdef __linux__
// linux code goes here.
#elif _WIN32
#include <winsock2.h>
#include <windows.h>
#include <psapi.h>
#include <iphlpapi.h>
#else
#error "OS not supported!"
#endif


#include "project.h"
#include "solver.h"

project::project() {
    task = NONE;
}

project::~project() {
}

void project::execute_task() {
    try {
        switch (task) {
        case LOAD_FES:
            std::cout << "Loading " << name << ".fes\n";
            model.frm.read_prj_file(full_path_name);
            model.sld.read_prj_file(full_path_name);
            model.msh.read_prj_file(full_path_name);
            break;
        case SAVE_FES:
            std::cout << "Saving " << name << ".fes\n";
            model.frm.write_prj_file(full_path_name);
            model.sld.write_prj_file(full_path_name);
            model.msh.write_prj_file(full_path_name);
            break;
        case LOAD_POLY:
            std::cout << "Loading " << name << ".poly\n";
            model.sld.read_poly_file(full_path_name);
            break;
        case SAVE_POLY:
            break;
        case LOAD_STL:
            std::cout << "Loading " << name << ".stl\n";
            model.sld.read_stl_file(full_path_name);
            model.create_tri_mesh();
            break;
        case LOAD_HFSS:
            std::cout << "Loading " << name << ".hfss\n";
            model.wrap_hfss(full_path_name, aux_path);
            break;
        case LOAD_AEDT:
            std::cout << "Loading " << name << ".aedt\n";
            model.wrap_aedt(full_path_name, aux_path);
            break;
        case RUN_TETGEN:
            std::cout << "Loading " << name << ".poly products\n";
            model.msh.read_tetgen_files(full_path_name);
            break;
        case RUN_TRIANGLE:
            std::cout << "Loading " << name << ".poly products\n";
            model.msh.read_triangle_files(full_path_name);
            break;
        case REFINE_HOMOGENEOUSLY:
            std::cout << "Refining homogeneously\n";
            model.msh.refine_homogeneous();
            break;
        case NONE:
            std::cout << "Nothing to do\n";
            break;
        }
        if (task == ANALYZE) {
            std::cout << "Running solver\n";
            solver sol(model);
        }
    } catch (std::string& str) {
        std::cout << "Error: " << str << "\n";
    }
}

std::string project::get_info() {
    std::stringstream tag;
    std::string name = get_var("COMPUTERNAME");
    if (name.size() != 0) {
        tag << "COMPUTERNAME         = " << name;
    } else {
        tag << "COMPUTERNAME         = ?";
    }
    tag << "\n";
    std::string cores = get_var("NUMBER_OF_PROCESSORS");
    if (cores.size() != 0) {
        tag << "NUMBER_OF_PROCESSORS = " << cores;
    } else {
        tag << "NUMBER_OF_PROCESSORS = ?";
    }
    tag << "\n";
    std::string threads = get_var("OMP_NUM_THREADS");
    if (threads.size() != 0) {
        tag << "OMP_NUM_THREADS      = " << threads;
    } else {
        tag << "OMP_NUM_THREADS      = " << cores; // automatically setting to the number of cores
    }
    tag << "\n";
    return "hello";//tag.str();
}

int project::get_num_proc() {
    return atoi(get_var("NUMBER_OF_PROCESSORS").data());
}

std::string project::get_var(const std::string name) {
    char* ptr = getenv(name.c_str());
    std::string ret;
    if (ptr == NULL) {
        ret = std::string("");
    } else {
        ret = std::string(ptr);
    }
    return ret;
}

int project::get_int(const std::string name) {
    const std::string data = get_var(name);
    int ret = -1;
    if (data.size() != 0) {
        ret = atoi(data.c_str());
    }
    return ret;
}

std::vector<std::string> project::get_mac() {
    std::vector<std::string> vMacAddresses;
    IP_ADAPTER_INFO AdapterInfo[32];        // Allocate information for up to 32 NICs
    DWORD dwBufLen = sizeof(AdapterInfo);   // Save memory size of buffer
    DWORD dwStatus = GetAdaptersInfo(       // Call GetAdapterInfo
                         AdapterInfo,       // [out] buffer to receive data
                         &dwBufLen);        // [in] size of receive data buffer
    // No network card? Other error?
    if (dwStatus != ERROR_SUCCESS) {
        return vMacAddresses;
    }
    PIP_ADAPTER_INFO pAdapterInfo = AdapterInfo;
    char szBuffer[512];
    while (pAdapterInfo) {
        if (pAdapterInfo->Type == MIB_IF_TYPE_ETHERNET) {
            snprintf(szBuffer, sizeof(szBuffer),
                     "%.2x-%.2x-%.2x-%.2x-%.2x-%.2x", pAdapterInfo->Address[0],
                     pAdapterInfo->Address[1], pAdapterInfo->Address[2],
                     pAdapterInfo->Address[3], pAdapterInfo->Address[4],
                     pAdapterInfo->Address[5]);
            vMacAddresses.push_back(std::string(szBuffer));
        }
        pAdapterInfo = pAdapterInfo->Next;
    }
    return vMacAddresses;
}

bool project::check_mac() {
    bool found = false;
    std::vector<std::string> vMacAddresses = get_mac();
    for (size_t i = 0; i < vMacAddresses.size(); i++) {
        if (vMacAddresses[i] == "00-00-00-00-00-00") {
            std::cout << "Welcome 00-00-00-00-00-00!\n";
            found |= true;
        }
    }
    return found;
}

std::string project::set_priority(unsigned int lvl) {
    HANDLE process = GetCurrentProcess();
    switch (lvl) {
    case 0:
        SetPriorityClass(process, NORMAL_PRIORITY_CLASS);
        break;
    case 1:
        SetPriorityClass(process, HIGH_PRIORITY_CLASS);
        break;
    case 2:
        SetPriorityClass(process, REALTIME_PRIORITY_CLASS);
        break;
    default:
        SetPriorityClass(process, NORMAL_PRIORITY_CLASS);
    }
    return get_priority();
}

std::string project::get_priority() {
    DWORD dwPriClass = GetPriorityClass(GetCurrentProcess());
    std::stringstream out;
    if (dwPriClass == REALTIME_PRIORITY_CLASS) {
        out << "REALTIME";
    } else if (dwPriClass == HIGH_PRIORITY_CLASS) {
        out << "HIGH";
    } else if (dwPriClass == NORMAL_PRIORITY_CLASS) {
        out << "NORMAL";
    } else {
        out << "level = " << dwPriClass;
    }
    out << " priority process";
    return out.str();
}

std::string project::get_proc_mem() {
    std::stringstream out;
    HANDLE hProcess = GetCurrentProcess();
    if (NULL == hProcess) {
        out << "Memory stats: Failed to acquire process handle";
        return out.str();
    }
    PROCESS_MEMORY_COUNTERS pmc;
    if (!GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc))) {
        out << "Memory stats: Failed to acquire process memory information";
    } else {
        double physiPeak = (double) pmc.PeakWorkingSetSize;
        double physiPres = (double) pmc.WorkingSetSize;
        out << "Used RAM: ";
        out << physiPres / 1048576;
        out << " MB |";
        out << physiPeak / 1048576;
        out << "|";
    }
    CloseHandle(hProcess);
    return out.str();
}

std::string project::get_sys_mem() {
    std::stringstream out;
#ifdef __linux__
    out << "Linux to be implemented" << std::endl;
#elif _WIN32
    MEMORYSTATUSEX statex;
    statex.dwLength = sizeof(statex);
    if (GlobalMemoryStatusEx(&statex)) {
        out << "Free RAM: ";
        out << statex.ullAvailPhys / 1048576;
        out << " MB";
    } else {
        out << "No memory information available";
    }
    out << "\n";
#else
#error "OS not supported!"
#endif
    return out.str();
}

std::string project::get_loc_time() {
    time_t ct = time(NULL);
    return std::string(asctime(localtime(&ct)));
}

std::string project::get_gmt() {
    time_t ct = time(NULL);
    return std::string(asctime(gmtime(&ct)));
}

std::string project::get_stats(timer& t) {
    return "- " /*+ get_proc_mem() + " - " */+ t.strtoc() + "\n";
}

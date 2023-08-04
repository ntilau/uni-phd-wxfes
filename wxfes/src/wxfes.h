#ifndef FES_H
#define FES_H

// #define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
// #define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);

#include <wx/wxprec.h>

#ifndef WX_PRECOMP
#include <wx/wx.h>
#endif

#include <wx/clipbrd.h>
#include <wx/cmdline.h>
#include <wx/laywin.h>
#include <wx/wfstream.h>
#include <wx/stdpaths.h>
#include <wx/string.h>
#include <wx/menu.h>
#include <wx/textctrl.h>
#include <wx/textdlg.h>
#include <wx/memory.h>
#include <wx/log.h>
#include <wx/treectrl.h>
#include <wx/filedlg.h>
#include <wx/choicdlg.h>
#include <wx/choice.h>
#include <wx/colordlg.h>
#include <wx/numdlg.h>
#include <wx/artprov.h>
#include <wx/image.h>
#include <wx/imaglist.h>
#include <wx/math.h>
#include <wx/renderer.h>
#include <wx/wupdlock.h>
#include <wx/filename.h>
#include <wx/thread.h>
#include <wx/progdlg.h>
#include <wx/event.h>
#include <wx/process.h>
#include <wx/dir.h>

#include "wxVTKRenderWindowInteractor.h"
#include <vtkVersion.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkActor2D.h>

// interactor styles
#include <vtkXYPlotWidget.h>
#include <vtkInteractorStyleDrawPolygon.h>
#include <vtkActor.h>
#include <vtkAffineRepresentation2D.h>
#include <vtkAffineWidget.h>
#include <vtkAppendPolyData.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkRenderView.h>
//
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkCoordinate.h>
#include <vtkAxesActor.h>
#include <vtkPointData.h>
#include <vtkLookupTable.h>
#include <vtkLine.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderView.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkTriangle.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include <future>         // std::async, std::future
#include <string>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef __linux__
#include <wx/wx.h>
#elif _WIN32
#include <windows.h>
#else
#error "OS not supported!"
#endif

#include <sstream>
#include <ostream>
#include <cstdlib>
#include "model.h"
#include "project.h"

wxDECLARE_EVENT(wxEVT_COMMAND_MYTHREAD_COMPLETED, wxThreadEvent);
wxDECLARE_EVENT(wxEVT_COMMAND_MYTHREAD_UPDATE, wxThreadEvent);

class MyApp: public wxApp {
public:
    virtual bool OnInit();
    /// Allows for command line entries to run solver to free GUI memory
    /// FES project file must be set up with proper model, analysis and required results
    virtual int OnExit();
    virtual int OnRun();
    virtual void OnInitCmdLine(wxCmdLineParser& parser);
    virtual bool OnCmdLineParsed(wxCmdLineParser& parser);
private:
    bool silent_mode;
};

static const wxCmdLineEntryDesc g_cmdLineDesc [] = {
    {
        wxCMD_LINE_SWITCH, "h", "help", "displays help on the command line parameters",
        wxCMD_LINE_VAL_NONE, wxCMD_LINE_OPTION_HELP
    },
//    {
//        wxCMD_LINE_SWITCH, "t", "test", "test switch",
//        wxCMD_LINE_VAL_NONE, wxCMD_LINE_PARAM_MANDATORY
//    },
    { wxCMD_LINE_OPTION, "s", "silent", "disables the GUI", wxCMD_LINE_VAL_STRING},
    { wxCMD_LINE_NONE }
};

DECLARE_APP(MyApp)

class MyFrame;

class MyThread : public wxThread {
public:
    MyThread(MyFrame *handler) : wxThread(wxTHREAD_DETACHED) {
        m_pHandler = handler;
    }
    ~MyThread();
protected:
    virtual ExitCode Entry();
    MyFrame *m_pHandler;
    timer t;
};

class material_dialog : public wxDialog {
public:
    material_dialog(wxWindow* parent, wxWindowID id, mdl_mtrl& mtrl);
    virtual ~material_dialog();
    void onCancel(wxCommandEvent& pEvent);
    void onOk(wxCommandEvent& pEvent);
private:
    wxTextCtrl* epsr, *mur, *sigma, *tand, *type, *name;
    mdl_mtrl &mtrl;
};

class boundary_dialog : public wxDialog {
public:
    boundary_dialog(wxWindow* parent, wxWindowID id, mdl_bc& bc, wxArrayString& str);
    virtual ~boundary_dialog();
    void onCancel(wxCommandEvent& pEvent);
    void onOk(wxCommandEvent& pEvent);
private:
    wxTextCtrl* name;
    wxChoice *type;
    mdl_bc& bc;
};


class MyFrame: public wxFrame {
public:
    MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
    virtual ~MyFrame();
    void CleanVTKFrame();
    void UpdateSystemMemory();
    void UpdateMesh();
    void UpdateSolid();
    void UpdateResults();
    void OnIdleEvent(wxIdleEvent& event);
    void OnQuit(wxCommandEvent& event);
    //void OnExit(wxCommandEvent& event);
    void OnAbout(wxCommandEvent& event);
    void OnOpen(wxCommandEvent& event);
    void OnSave(wxCommandEvent& event);
    void OnImportPoly(wxCommandEvent& event);
    void OnImportHfss(wxCommandEvent& event);
    void OnImportAedt(wxCommandEvent& event);
    void run_analysis(wxCommandEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnSashDrag(wxSashEvent& event);
    void OnClearLog(wxCommandEvent& event);
    void OnSetBackgroundColor(wxCommandEvent& event);
    // Tree_Ctrl
    void OnBeginDrag(wxTreeEvent& event);
    void OnBeginRDrag(wxTreeEvent& event);
    void OnEndDrag(wxTreeEvent& event);
    void OnBeginLabelEdit(wxTreeEvent& event);
    void OnEndLabelEdit(wxTreeEvent& event);
    void OnDeleteItem(wxTreeEvent& event);
    void OnGetInfo(wxTreeEvent& event);
    void OnSetInfo(wxTreeEvent& event);
    void OnItemExpanded(wxTreeEvent& event);
    void OnItemExpanding(wxTreeEvent& event);
    void OnItemCollapsed(wxTreeEvent& event);
    void OnItemCollapsing(wxTreeEvent& event);
    void OnSelChanged(wxTreeEvent& event);
    void OnSelChanging(wxTreeEvent& event);
    void OnTreeKeyDown(wxTreeEvent& event);
    void OnItemActivated(wxTreeEvent& event);
    void OnItemStateClick(wxTreeEvent& event);
    void OnItemMenu(wxTreeEvent& event);
    void OnItemRClick(wxTreeEvent& event);
    // viewer
    void view_sld();
    void view_msh();
    // triangle & tetgen
    void setup_tetgen(wxCommandEvent& event);
    void setup_triangle(wxCommandEvent& event);
    void refine_homogeneously(wxCommandEvent& event);
    void set_model();

    // threads
    void DoStartThread();
    void DoPauseThread();
    void OnThreadUpdate(wxCommandEvent&);
    void OnThreadCompletion(wxCommandEvent&);
    void OnClose(wxCloseEvent&);
protected:
    MyThread *m_pThread = NULL;
    wxCriticalSection m_pThreadCS;    // protects the m_pThread pointer
    friend class MyThread;            // allow it to access our m_pThread
    wxString CurrentDocPath;
    wxTextCtrl *m_logCtrl;
    wxStreamToTextRedirector* redirect;
private:
    wxSashLayoutWindow *m_pSashWindowLeft;
    wxSashLayoutWindow *m_pSashWindowRight;
    wxSashLayoutWindow *m_pSashWindowBottom;
///// threading
//    wxThread *m_pThread;
/// vtk
    wxVTKRenderWindowInteractor* m_pVTKWindow;
    vtkRenderer* pRenderer;
    vtkRenderWindow* pRenderWindow;
    vtkSmartPointer<vtkActor> pSldActor, pMshActor, pFldActor;
    vtkSmartPointer<vtkScalarBarActor> scalarBar;
    vtkSmartPointer<vtkAxesActor> pAxesActor;
    double current_visibility = 1.0;
    std::vector<int> bcs_to_view;
    unsigned int bcs_default_size;
/// tree
    wxTreeCtrl *m_treeCtrl;
    wxTreeItemId tree_project,
                 tree_model,
                 tree_manifold, tree_boundaries, tree_materials,
                 tree_mesh,
                 tree_analysis,
                 tree_results;
    std::map<wxItemId<void*>::Type, unsigned int> bc_setting, mtrl_setting;
    wxItemId<void*>::Type prj_id;
    project prj;
    DECLARE_EVENT_TABLE()
};

/////////////////////////////
/// STARTING OPERARATIONS ///
/////////////////////////////
#define MY_FRAME      101
#define MY_VTK_WINDOW 102
#define MY_WINDOW_LEFT  103
#define MY_WINDOW_RIGHT 104
#define MY_WINDOW_BOTTOM 105
#define MY_TREE_CTRL 110

enum {
    FES_Quit = wxID_EXIT,
    FES_About = wxID_ABOUT,
    FES_Open = wxID_OPEN,
    FES_Save = wxID_SAVE,
    FES_Tree_Ctrl = MY_TREE_CTRL,
    FES_SetBGColor,
    FES_Import_poly,
    FES_Import_hfss,
    FES_Import_aedt,
    FES_Run,
    FES_Tetgen,
    FES_Triangle,
    FES_RefineHomogeneously,
};

#endif // FES_H

#include "wxfes.h"
#include <vector>
#include <list>

bool MyApp::OnInit() {
//    if (!wxApp::OnInit())
//        return false;
    MyFrame *frame = new MyFrame(_T("Finite Element Software"), wxDefaultPosition, wxSize(800, 600));
    //frame->SetIcon(wxICON(IDR_MAINFRAME));
    frame->Center();
    frame->Show(true);
    frame->Maximize();
    //SetTopWindow(frame);
    return true;
}

int MyApp::OnExit() {
    // clean up
    return 0;
}

int MyApp::OnRun() {
    int exitcode = wxApp::OnRun();
    wxTheClipboard->Flush();
    if (exitcode!=0)
        return exitcode;
}

void MyApp::OnInitCmdLine(wxCmdLineParser& parser) {
    parser.SetDesc(g_cmdLineDesc);
    // must refuse '/' as parameter starter or cannot use "/path" style paths
    parser.SetSwitchChars(wxT("-"));
}

bool MyApp::OnCmdLineParsed(wxCmdLineParser& parser) {
    silent_mode = parser.Found(wxT("s"));
    // to get at your unnamed parameters use
    wxArrayString files;
    for (int i = 0; i < parser.GetParamCount(); i++) {
        files.Add(parser.GetParam(i));
    }
    // and other command line parameters

    // then do what you need with them.

    /// run the task -- cmd prompt missing...

    return true;
}

material_dialog::material_dialog(wxWindow* parent, wxWindowID id, mdl_mtrl& _mtrl)
    : wxDialog(parent, id, "Set material", wxDefaultPosition, wxDefaultSize, wxDEFAULT_DIALOG_STYLE),
      mtrl(_mtrl) {
    wxPanel* panel = new wxPanel(this, wxID_ANY);
    epsr = new wxTextCtrl(panel, wxID_ANY, wxString::Format(wxT("%e"), mtrl.epsr));
    mur = new wxTextCtrl(panel, wxID_ANY, wxString::Format(wxT("%e"), mtrl.mur));
    sigma = new wxTextCtrl(panel, wxID_ANY, wxString::Format(wxT("%e"), mtrl.sigma));
    tand = new wxTextCtrl(panel, wxID_ANY, wxString::Format(wxT("%e"), mtrl.tand));
    type = new wxTextCtrl(panel, wxID_ANY, mtrl.type);
    name = new wxTextCtrl(panel, wxID_ANY, mtrl.name);
    wxButton* cancelButton = new wxButton(panel, wxID_ANY, "Cancel");
    wxButton* okButton = new wxButton(panel, wxID_ANY, "Ok");
    wxGridSizer *sizer = new wxGridSizer(7, 2, 3, 3);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "Solid:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(name, 1, wxALL);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "Type:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(type, 1, wxALL);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "EpsR:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(epsr, 1, wxALL);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "MuR:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(mur, 1, wxALL);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "Sigma:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(sigma, 1, wxALL);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "tanDelta:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(tand, 1, wxALL);
    sizer->Add(okButton, 1, wxALIGN_CENTER_HORIZONTAL|wxALIGN_BOTTOM);
    sizer->Add(cancelButton, 1, wxALIGN_CENTER_HORIZONTAL|wxALIGN_BOTTOM);
    panel->SetSizer(sizer);
    Bind(wxEVT_COMMAND_BUTTON_CLICKED, &material_dialog::onCancel, this, cancelButton->GetId());
    Bind(wxEVT_COMMAND_BUTTON_CLICKED, &material_dialog::onOk, this, okButton->GetId());
    Center();
}

material_dialog::~material_dialog() {}

void material_dialog::onCancel(wxCommandEvent& WXUNUSED(pEvent))    {
    EndModal(wxID_CANCEL);
    Destroy();
}

void material_dialog::onOk(wxCommandEvent& WXUNUSED(pEvent))    {
    mtrl.name = name->GetValue().ToStdString();
    mtrl.type = type->GetValue().ToStdString();
    mtrl.epsr = atof(epsr->GetValue().ToStdString().c_str());
    mtrl.mur = atof(mur->GetValue().ToStdString().c_str());
    mtrl.sigma = atof(sigma->GetValue().ToStdString().c_str());
    mtrl.tand = atof(tand->GetValue().ToStdString().c_str());
    EndModal(wxID_OK);
    Destroy();
}

boundary_dialog::boundary_dialog(wxWindow* parent, wxWindowID id, mdl_bc& _bc, wxArrayString& str)
    : wxDialog(parent, id, "Set boundary condition", wxDefaultPosition, wxSize(250,120), wxDEFAULT_DIALOG_STYLE),
      bc(_bc) {
    wxPanel* panel = new wxPanel(this, wxID_ANY);
    name = new wxTextCtrl(panel, wxID_ANY, bc.name);
    type = new wxChoice(panel, wxID_ANY, wxDefaultPosition, wxDefaultSize, str);
    for(unsigned int i = 0; i< str.GetCount(); i++)
        if(strcmp(bc.type.data(),str[i].c_str()) == 0)
            type->SetSelection(i);
    wxButton* cancelButton = new wxButton(panel, wxID_ANY, "Cancel");
    wxButton* okButton = new wxButton(panel, wxID_ANY, "Ok");
    wxGridSizer *sizer = new wxGridSizer(3, 2, 3, 3);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "Name:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(name, 1, wxALL);
    sizer->Add(new wxStaticText(panel, wxID_ANY, "Type:"), 1, wxALIGN_CENTER_HORIZONTAL);
    sizer->Add(type, 1, wxALL);
    sizer->Add(okButton, 1, wxALIGN_CENTER_HORIZONTAL|wxALIGN_BOTTOM);
    sizer->Add(cancelButton, 1, wxALIGN_CENTER_HORIZONTAL|wxALIGN_BOTTOM);
    panel->SetSizer(sizer);
    Bind(wxEVT_COMMAND_BUTTON_CLICKED, &boundary_dialog::onCancel, this, cancelButton->GetId());
    Bind(wxEVT_COMMAND_BUTTON_CLICKED, &boundary_dialog::onOk, this, okButton->GetId());
    Center();
}

boundary_dialog::~boundary_dialog() {}

void boundary_dialog::onCancel(wxCommandEvent& WXUNUSED(pEvent)) {
    EndModal(wxID_CANCEL);
    Destroy();
}

void boundary_dialog::onOk(wxCommandEvent& WXUNUSED(pEvent)) {
    bc.name = name->GetValue().ToStdString();
    bc.type = type->GetStringSelection().ToStdString();
    if(strcmp(type->GetStringSelection().ToStdString().data(), "Voltage") == 0) {
        std::stringstream ss;
        ss << bc.voltage;
        wxTextEntryDialog dlg(nullptr, "Enter voltage [V]", "Boundary parameter", ss.str());
        if(dlg.ShowModal() == wxID_OK) {
            bc.voltage = atof(dlg.GetValue().ToStdString().c_str());
        }
    }
    EndModal(wxID_OK);
    Destroy();
}

MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size) :
    wxFrame((wxFrame *) NULL, -1, title, pos, size), m_pSashWindowLeft(0), m_pSashWindowRight(0),
    m_pSashWindowBottom(0), m_treeCtrl(NULL), pSldActor(NULL) {
    // DisableProcessWindowsGhosting();
    wxMenu* menuFile = new wxMenu(_T(""), wxMENU_TEAROFF);
    wxMenu* menuView = new wxMenu(_T(""), wxMENU_TEAROFF);
    wxMenu* menuRun = new wxMenu(_T(""), wxMENU_TEAROFF);
    wxMenu* helpMenu = new wxMenu;
    helpMenu->Append(FES_About, _T("&About...\tCtrl-A"), _T("Show about dialog"));
    menuFile->Append(FES_Open, _T("&Open\tAlt-O"), _T("Open project file"));
    menuFile->Append(FES_Save, _T("&Save\tAlt-S"), _T("Save project file"));
    menuFile->AppendSeparator();
    menuFile->Append(FES_Quit, _T("E&xit\tAlt-X"), _T("Quit FES"));
    menuView->Append(FES_SetBGColor, _T("&Background Color\tAlt-B"), _T("Set background color"));
    wxMenuBar* menuBar = new wxMenuBar();
    menuBar->Append(menuFile, _T("&File"));
    menuBar->Append(menuView, _T("&View"));
    menuBar->Append(helpMenu, _T("&Help"));
    SetMenuBar(menuBar);
    CreateStatusBar(3);
    SetStatusText(_T("Welcome to FES"),0);
    UpdateSystemMemory();
    m_pSashWindowBottom = new wxSashLayoutWindow(this, MY_WINDOW_BOTTOM,
            wxDefaultPosition, wxSize(size.GetWidth(), size.GetHeight()*.2),
            wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    m_pSashWindowBottom->SetDefaultSize(wxSize(size.GetWidth(),size.GetHeight()*.2));
    m_pSashWindowBottom->SetOrientation(wxLAYOUT_HORIZONTAL);
    m_pSashWindowBottom->SetAlignment(wxLAYOUT_BOTTOM);
    m_pSashWindowBottom->SetBackgroundColour(*wxWHITE);
    m_pSashWindowBottom->SetSashVisible(wxSASH_TOP, TRUE);
    m_pSashWindowLeft = new wxSashLayoutWindow(this, MY_WINDOW_LEFT,
            wxDefaultPosition, wxSize(size.GetWidth()*.25, size.GetHeight()*.8),
            wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    m_pSashWindowLeft->SetDefaultSize(wxSize(size.GetWidth()*.25,size.GetHeight()*.8));
    m_pSashWindowLeft->SetOrientation(wxLAYOUT_VERTICAL);
    m_pSashWindowLeft->SetAlignment(wxLAYOUT_LEFT);
    m_pSashWindowLeft->SetSashVisible(wxSASH_RIGHT, TRUE);
    m_pSashWindowRight = new wxSashLayoutWindow(this, MY_WINDOW_RIGHT,
            wxDefaultPosition, wxSize(size.GetWidth()*.75, size.GetHeight()*.8),
            wxNO_BORDER | wxSW_3D | wxCLIP_CHILDREN);
    m_pSashWindowRight->SetDefaultSize(wxSize(size.GetWidth()*.75,size.GetHeight()*.8));
    m_pSashWindowRight->SetOrientation(wxLAYOUT_VERTICAL);
    m_pSashWindowRight->SetAlignment(wxLAYOUT_RIGHT);

    m_logCtrl = new wxTextCtrl(m_pSashWindowBottom, wxID_ANY, wxEmptyString,
                               wxDefaultPosition,
                               wxSize(size.GetWidth(),size.GetHeight()*.2),
                               wxTE_MULTILINE|wxTE_PROCESS_ENTER|wxTE_READONLY,
                               wxDefaultValidator, "Logger");
    redirect = new wxStreamToTextRedirector(m_logCtrl);
#ifdef _WIN32
    m_pVTKWindow = new wxVTKRenderWindowInteractor(m_pSashWindowRight, MY_VTK_WINDOW);
    pRenderer = vtkRenderer::New();
    pRenderer->GetActiveCamera()->ParallelProjectionOn();
    m_pVTKWindow->GetRenderWindow()->AddRenderer(pRenderer);
    //m_pVTKWindow->SetInteractorStyle(vtkSmartPointer<vtkInteractorStyleDrawPolygon>::New());
    vtkInteractorStyleSwitch::SafeDownCast(m_pVTKWindow->GetInteractorStyle())->AutoAdjustCameraClippingRangeOn();
    vtkInteractorStyleSwitch::SafeDownCast(m_pVTKWindow->GetInteractorStyle())->SetCurrentStyleToTrackballCamera();
//    pRenderer->SetBackground(double(105.0/255.0),
//                             double(144.0/255.0),
//                             double(198.0/255.0));
    pRenderer->SetBackground(double(240.0/255.0),
                             double(240.0/255.0),
                             double(240.0/255.0));
    /*
    vtkSmartPointer<vtkRenderView> renderView = vtkSmartPointer<vtkRenderView>::New();
    renderView->SetInteractionMode(vtkRenderView::INTERACTION_MODE_3D);
    renderView->SetInteractor(m_pVTKWindow);
    renderView->Update();*/
    pAxesActor = vtkSmartPointer<vtkAxesActor>::New();
    pRenderer->AddActor(pAxesActor);
    pRenderer->ResetCamera();
    m_pVTKWindow->Render();
#endif

    CurrentDocPath = wxStandardPaths::Get().GetLocalDataDir();
    m_treeCtrl = new wxTreeCtrl(m_pSashWindowLeft, MY_TREE_CTRL, wxDefaultPosition,
                                wxDefaultSize, wxTR_DEFAULT_STYLE);
    m_treeCtrl->AddRoot("Project");
    tree_project = m_treeCtrl->GetRootItem();
    m_treeCtrl->AppendItem(tree_project, "Model");
    tree_model = m_treeCtrl->GetLastChild(tree_project);
    m_treeCtrl->AppendItem(tree_model, "Boundaries");
    tree_boundaries = m_treeCtrl->GetLastChild(tree_model);
    m_treeCtrl->AppendItem(tree_model, "Materials");
    tree_materials = m_treeCtrl->GetLastChild(tree_model);
    m_treeCtrl->AppendItem(tree_project, "Mesh");
    tree_mesh = m_treeCtrl->GetLastChild(tree_project);
    m_treeCtrl->AppendItem(tree_project, "Analysis");
    tree_analysis = m_treeCtrl->GetLastChild(tree_project);
    m_treeCtrl->AppendItem(tree_project, "Results");
    tree_results = m_treeCtrl->GetLastChild(tree_project);
    m_treeCtrl->ExpandAllChildren(tree_project);
    // ADDING FILE TO THE PATH in order to manage tetgen & triangle
    // wxString set_path("setx PATH ");
    // set_path.append(wxStandardPaths::Get().GetExecutablePath().BeforeLast(wxFILE_SEP_PATH));
    // wxExecute(set_path, wxEXEC_HIDE_CONSOLE);
}

MyFrame::~MyFrame() {
    if (m_pVTKWindow)
        m_pVTKWindow->Delete();
    if (pRenderer != NULL)
        pRenderer->Delete();
    if (m_treeCtrl != NULL)
        delete m_treeCtrl;
    if (m_logCtrl != NULL)
        delete m_logCtrl;
    if (redirect != NULL)
        delete redirect;
}

void MyFrame::UpdateSystemMemory() {
    SetStatusText(project().get_proc_mem(),1);
    SetStatusText(project().get_sys_mem(),2);
}

/// event handlers
void MyFrame::OnIdleEvent(wxIdleEvent& event) {
//    std::cout << "idle\n";
//    UpdateMesh();
}


void MyFrame::OnQuit(wxCommandEvent& WXUNUSED(event)) {
    Close(TRUE);
}

void MyFrame::OnAbout(wxCommandEvent& WXUNUSED(event)) {
    wxString msg;
    msg.Printf(_T("wxWidgets+VTK experiment provided\nby Laurent Ntibarikure since 2015"));
    wxMessageBox(msg, _T("About Finite Element Software"), wxOK | wxICON_INFORMATION, this);
}

void MyFrame::OnSetBackgroundColor(wxCommandEvent& WXUNUSED(event)) {
    wxColourData data;
    data.SetChooseFull(true);
    for(int i = 0; i < 16; i++) {
        wxColour colour(i*16, i*16, i*16);
        data.SetCustomColour(i, colour);
    }
    wxColourDialog dlg(this, &data);
    dlg.SetTitle("Viewer background color");
    if(dlg.ShowModal() == wxID_OK) {
        wxColourData retData = dlg.GetColourData();
        wxColour col = retData.GetColour();
        wxInt32 RGB = col.GetRGB();
        int R = (0x000000FF & RGB);
        int G = (0x0000FF00 & RGB) >> 8;
        int B = (0x00FF0000 & RGB) >> 16;
        pRenderer->SetBackground(double(R)/255.0,
                                 double(G)/255.0,
                                 double(B)/255.0);
        m_pVTKWindow->Render();
    }
}


void MyFrame::OnOpen(wxCommandEvent& WXUNUSED(event)) {
    wxFileDialog* openFileDialog = new wxFileDialog(this, _("Open FES file"), "", "",
            "FES files (*.fes)|*.fes", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
    if (openFileDialog->ShowModal() == wxID_OK) {
        CurrentDocPath = openFileDialog->GetDirectory();
        wxFileName name(openFileDialog->GetFilename());
        prj.name = name.GetName();
        prj.data_path = openFileDialog->GetDirectory().Append(wxFILE_SEP_PATH);
        prj.full_path_name = prj.data_path + prj.name;
        m_treeCtrl->SetItemText(tree_project, std::string("Project: \"" + prj.name + "\"").data());
        prj.task = project::LOAD_FES;
        DoStartThread();
    }
    openFileDialog->Destroy();
}

void MyFrame::OnSave(wxCommandEvent& WXUNUSED(event)) {
    wxFileDialog* saveFileDialog = new wxFileDialog(this, _("Save FES file"), "", "",
            "FES files (*.fes)|*.fes", wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
    saveFileDialog->SetFilename(std::string(prj.name + ".fes").data() );
    if (saveFileDialog->ShowModal() == wxID_OK) {
        CurrentDocPath = saveFileDialog->GetDirectory();
        wxFileName name(saveFileDialog->GetFilename());
        prj.name = name.GetName();
        prj.data_path = saveFileDialog->GetDirectory().Append(wxFILE_SEP_PATH);
        prj.full_path_name = prj.data_path + prj.name;
        m_treeCtrl->SetItemText(tree_project, std::string("Project: \"" + prj.name + "\"").data());
        prj.task = project::SAVE_FES;
        DoStartThread();
    }
    saveFileDialog->Destroy();
}

void MyFrame::OnSize(wxSizeEvent& event) {
    if (m_pSashWindowLeft != 0)
        m_pSashWindowLeft->SetDefaultSize(
            wxSize(event.GetSize().GetWidth() / 4,
                   event.GetSize().GetHeight()));
    if (m_pSashWindowRight != 0)
        m_pSashWindowRight->SetDefaultSize(
            wxSize(event.GetSize().GetWidth() / 4 * 3,
                   event.GetSize().GetHeight()));
    wxLayoutAlgorithm layout;
    layout.LayoutFrame(this);
}

void MyFrame::OnSashDrag(wxSashEvent& event) {
    if (event.GetDragStatus() == wxSASH_STATUS_OUT_OF_RANGE)
        return;
    switch (event.GetId()) {
    case MY_WINDOW_BOTTOM: {
        if (m_pSashWindowBottom != 0)
            m_pSashWindowBottom->SetDefaultSize(
                wxSize(GetSize().GetWidth(), event.GetDragRect().height));
        break;
    }
    case MY_WINDOW_LEFT: {
        if (m_pSashWindowLeft != 0)
            m_pSashWindowLeft->SetDefaultSize(
                wxSize(event.GetDragRect().width, GetSize().GetHeight()));
        if (m_pSashWindowRight != 0)
            m_pSashWindowRight->SetDefaultSize(
                wxSize(GetSize().GetWidth() - event.GetDragRect().width,
                       GetSize().GetHeight()));
        break;
    }
    case MY_WINDOW_RIGHT: {
        if (m_pSashWindowRight != 0)
            m_pSashWindowRight->SetDefaultSize(
                wxSize(event.GetDragRect().width, GetSize().GetHeight()));
        if (m_pSashWindowLeft != 0)
            m_pSashWindowLeft->SetDefaultSize(
                wxSize(GetSize().GetWidth() - event.GetDragRect().width,
                       GetSize().GetHeight()));
        break;
    }
    }
    wxLayoutAlgorithm layout;
    layout.LayoutFrame(this);
}

void MyFrame::OnItemRClick(wxTreeEvent& event) {
    wxTreeItemId item = event.GetItem();
    if(m_treeCtrl->GetItemParent(item).IsOk()) {
        if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Model") == 0) {
            wxMenu menu;
            menu.Append(FES_Import_poly, wxT("Import *.poly"));
            menu.Append(FES_Import_hfss, wxT("Import *.hfss"));
            menu.Append(FES_Import_aedt, wxT("Import *.aedt"));
            PopupMenu(&menu, event.GetPoint());
        }
        if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Mesh") == 0) {
            wxMenu menu;
            menu.Append(FES_Tetgen, wxT("Run Tetgen"));
            menu.Append(FES_Triangle, wxT("Run Triangle"));
            menu.Append(FES_RefineHomogeneously, wxT("Refine homogeneously"));
            PopupMenu(&menu, event.GetPoint());
        }
        if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Analysis") == 0) {
            wxMenu menu;
            menu.Append(FES_Run, wxT("Run"));
            PopupMenu(&menu, event.GetPoint());
        }
        if(strcmp(m_treeCtrl->GetItemText(m_treeCtrl->GetItemParent(item)).c_str(), "Boundaries") == 0) {
            unsigned int frm_id = 0;
            for(unsigned int i=0; i< prj.model.frm.frm_type.size(); i++)
                for(std::list<std::string>::iterator it = prj.model.frm.frm_type[i].begin(); it != prj.model.frm.frm_type[i].end(); it++)
                    if(strcmp(prj.model.frm.type.data(), it->data()) == 0)
                        frm_id = i;
            wxArrayString str;
            for(std::list<std::string>::iterator it = prj.model.frm.bc_type[frm_id].begin(); it != prj.model.frm.bc_type[frm_id].end(); it++)
                str.Add(wxString(*it));
            boundary_dialog bc_dlg(nullptr, wxID_ANY, prj.model.frm.bcs[bc_setting[item.GetID()]], str);
            if(bc_dlg.ShowModal() == wxID_OK) {
                std::stringstream text;

                text << prj.model.frm.bcs[bc_setting[item.GetID()]].label << " - " <<
                     prj.model.frm.bcs[bc_setting[item.GetID()]].name << " - " <<
                     prj.model.frm.bcs[bc_setting[item.GetID()]].type;
                m_treeCtrl->SetItemText(item,text.str());
            }
        }
        if(strcmp(m_treeCtrl->GetItemText(m_treeCtrl->GetItemParent(item)).c_str(), "Materials") == 0) {
            material_dialog mtrl_dlg(nullptr, wxID_ANY, prj.model.frm.mtrls[mtrl_setting[item.GetID()]]);
            if(mtrl_dlg.ShowModal() == wxID_OK) {
                std::stringstream text;
                text << prj.model.frm.mtrls[mtrl_setting[item.GetID()]].label << " - " <<
                     prj.model.frm.mtrls[mtrl_setting[item.GetID()]].name << " - " <<
                     prj.model.frm.mtrls[mtrl_setting[item.GetID()]].type;
                m_treeCtrl->SetItemText(item,text.str());
            }
        }
        if(strcmp(m_treeCtrl->GetItemText(m_treeCtrl->GetItemParent(item)).c_str(), "Analysis") == 0) {
            wxArrayString str;
            for(unsigned int i=0; i< prj.model.frm.frm_type.size(); i++)
                for(std::list<std::string>::iterator it = prj.model.frm.frm_type[i].begin(); it != prj.model.frm.frm_type[i].end(); it++)
                    str.Add(wxString(*it));
            wxSingleChoiceDialog* bc_dlg = new wxSingleChoiceDialog(this, "Select one entry", "Set formulation", str);
            if (bc_dlg->ShowModal() == wxID_OK) {
                prj.model.frm.type = str[bc_dlg->GetSelection()];
                m_treeCtrl->SetItemText(item,prj.model.frm.type);
            }
        }
    } else {
        wxTextEntryDialog prj_dlg(this, "Set project name", "Project parameters", prj.name.data());
        if (prj_dlg.ShowModal() == wxID_OK) {
            prj.name = prj_dlg.GetValue().ToStdString();
        }
        if(prj.name.size() > 0)
            m_treeCtrl->SetItemText(tree_project, std::string("Project: \"" + prj.name + "\"").data());
    }
    UpdateSystemMemory();
}

void MyFrame::set_model() {
    if(prj.model.frm.bcs.size() == 0) {
        prj.model.frm.bcs.resize(prj.model.sld.bc_markers.size());
        for(unsigned int i=0; i<prj.model.sld.bc_markers.size(); i++)
            prj.model.frm.bcs[i].label = prj.model.sld.bc_markers[i];
    }
    if(prj.model.frm.mtrls.size() == 0) {
        prj.model.frm.mtrls.resize(prj.model.sld.mtrl_markers.size());
        for(unsigned int i=0; i<prj.model.sld.mtrl_markers.size(); i++)
            prj.model.frm.mtrls[i].label = prj.model.sld.mtrl_markers[i];
    }
    m_treeCtrl->DeleteChildren(tree_boundaries);
    bcs_to_view = prj.model.msh.fac_lab;
    bcs_to_view.erase(unique(bcs_to_view.begin(),bcs_to_view.end()), bcs_to_view.end());
    bcs_default_size = bcs_to_view.size();
    for(unsigned int i=0; i < prj.model.frm.bcs.size(); i++) {
        std::stringstream labs;
        labs << prj.model.frm.bcs[i].label <<
             " - " << prj.model.frm.bcs[i].name <<
             " - " << prj.model.frm.bcs[i].type;
        m_treeCtrl->AppendItem(tree_boundaries, labs.str());
        bc_setting[m_treeCtrl->GetLastChild(tree_boundaries).GetID()] = i;
    }
    m_treeCtrl->DeleteChildren(tree_materials);
    for(unsigned int i=0; i<prj.model.frm.mtrls.size(); i++) {
        std::stringstream labs;
        labs << prj.model.frm.mtrls[i].label <<
             " - " << prj.model.frm.mtrls[i].name <<
             " - " << prj.model.frm.mtrls[i].type;
        m_treeCtrl->AppendItem(tree_materials, labs.str());
        mtrl_setting[m_treeCtrl->GetLastChild(tree_materials).GetID()] = i;
    }
    m_treeCtrl->DeleteChildren(tree_analysis);
    m_treeCtrl->AppendItem(tree_analysis, prj.model.frm.type);
    m_treeCtrl->ExpandAllChildren(tree_project);
}

void MyFrame::OnItemActivated(wxTreeEvent& event) {
    wxTreeItemId item = event.GetItem();
    if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Model") == 0) {
        m_treeCtrl->ExpandAllChildren(item);
        if(current_visibility == 1.0)
            current_visibility = 0.5;
        else
            current_visibility = 1.0;
        set_model();
        view_sld();
    }
    if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Boundaries") == 0) {
        m_treeCtrl->ExpandAllChildren(item);
    }
    for(unsigned int i = 0; i < prj.model.frm.bcs.size(); i++) {
        //if(bc_setting[item.GetID()])
        std::stringstream labs;
        labs << prj.model.frm.bcs[i].label <<
             " - " << prj.model.frm.bcs[i].name <<
             " - " << prj.model.frm.bcs[i].type;
        if(strcmp(m_treeCtrl->GetItemText(item).c_str(), labs.str().c_str()) == 0) {
//            if(bcs_to_view.size() == bcs_default_size)
//                bcs_to_view.clear();
            if(std::find(bcs_to_view.begin(), bcs_to_view.end(), prj.model.frm.bcs[i].label) != bcs_to_view.end())
                bcs_to_view.erase(std::remove(bcs_to_view.begin(), bcs_to_view.end(), prj.model.frm.bcs[i].label), bcs_to_view.end());
            else
                bcs_to_view.push_back(prj.model.frm.bcs[i].label);
            view_sld();
        }
    }
    if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Mesh") == 0) {
        m_treeCtrl->ExpandAllChildren(item);
        if(current_visibility == 1.0)
            current_visibility = 0.5;
        else
            current_visibility = 1.0;
        set_model();
        //view_sld();
        view_msh();
    }
    if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Analysis") == 0) {
        m_treeCtrl->ExpandAllChildren(item);
    }
    if(strcmp(m_treeCtrl->GetItemText(item).c_str(), "Results") == 0) {
        m_treeCtrl->ExpandAllChildren(item);
        CleanVTKFrame();
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(size_t i = 0; i< prj.model.msh.nod_pos.size(); i++)
            points->InsertNextPoint(prj.model.msh.nod_pos[i][0], prj.model.msh.nod_pos[i][1], prj.model.msh.nod_pos[i][2]);
        vtkSmartPointer<vtkCellArray> faces =  vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkDoubleArray> faces_val =  vtkSmartPointer<vtkDoubleArray>::New();
        faces_val->SetNumberOfComponents(1);
        faces_val->SetName("Potential");
        for(size_t i = 0; i < prj.model.msh.fac_nodes.size(); i++) {
            vtkSmartPointer<vtkTriangle> face = vtkSmartPointer<vtkTriangle>::New();
            for(size_t j = 0; j < prj.model.msh.fac_nodes[i].size(); j++) {
                face->GetPointIds()->SetId(j, prj.model.msh.fac_nodes[i][j]);
            }
            faces->InsertNextCell(face);
//                faces_val->InsertNextValue(double(prj.model.msh.edges_marker[i][0])/13);
        }
        double max_val = -DBL_MAX, min_val = DBL_MAX;
        for(size_t i = 0; i < prj.model.msh.nod_pos.size(); i++) {
            max_val = std::max(max_val,prj.model.frm.sol_real[0][i]);
            min_val = std::min(min_val,prj.model.frm.sol_real[0][i]);
        }
        for(size_t i = 0; i < prj.model.msh.nod_pos.size(); i++) {
            faces_val->InsertNextValue((prj.model.frm.sol_real[0][i] - min_val) /(max_val - min_val));
        }
        vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
        polygonPolyData->SetPoints(points);
        polygonPolyData->SetPolys(faces);
        polygonPolyData->GetPointData()->SetScalars(faces_val);
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polygonPolyData);
        mapper->ScalarVisibilityOn();
        mapper->SetScalarModeToUsePointData();
        mapper->SetColorModeToMapScalars();

        scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar->SetLookupTable(mapper->GetLookupTable());
        scalarBar->SetTitle("Potential");
        scalarBar->SetNumberOfLabels(6);
        scalarBar->SetHeight(0.10);
        scalarBar->SetWidth(0.90);
        scalarBar->SetPosition(0.05, 0.05);
        scalarBar->SetOrientationToHorizontal();
        scalarBar->SetLabelFormat("%.2g");
        scalarBar->VisibilityOn();
        scalarBar->GetLabelTextProperty()->SetFontSize(6);

        // Create a lookup table to share between the mapper and the scalarbar
        vtkSmartPointer<vtkLookupTable> lut =
            vtkSmartPointer<vtkLookupTable>::New();
        lut->SetTableRange(min_val, max_val);
        lut->SetHueRange(0.667,0.0);
//        lut->SetSaturationRange (.5, .5);
//        lut->SetValueRange(1.0, 1.0);
//        lut->SetTableRange(0, 280);
        lut->SetNumberOfColors(256);

        lut->Build();

        mapper->SetLookupTable( lut );
        scalarBar->SetLookupTable( lut );
//            mapper->SetScalarModeToUseCellFieldData();
//            mapper->SetLookupTable(lut);
        pFldActor = vtkSmartPointer<vtkActor>::New();
        pFldActor->SetMapper(mapper);
//            pFldActor->GetProperty()->SetRepresentationToWireframe();
        pRenderer->AddActor(pFldActor);
        pRenderer->AddActor2D(scalarBar);
        double geom_dim = prj.model.msh.get_geom_dim();
        pAxesActor->SetTotalLength(geom_dim,geom_dim,geom_dim);
        prj.model.msh.get_bounding_info();
        std::vector<std::vector<double> > bbox = prj.model.msh.bounding_box;
        pRenderer->ResetCamera(bbox[0][0], bbox[1][0],
                               bbox[0][1], bbox[1][1],
                               bbox[0][2], bbox[1][2]);
        m_pVTKWindow->Render();
    }

    UpdateSystemMemory();
}

void MyFrame::setup_tetgen(wxCommandEvent& WXUNUSED(event)) {
    wxTextEntryDialog* setupTetgenDialog = new wxTextEntryDialog(this,
            "Command line options:\npq__a__AriYMS__T__dzjo_fengGOJBNEFICQVvh", "Setup Tetgen");
    setupTetgenDialog->SetValue (prj.model.sld.tetgen_command);
    if (prj.model.sld.dim == 3 && setupTetgenDialog->ShowModal() == wxID_OK) {
        prj.model.sld.tetgen_command = setupTetgenDialog->GetValue();
        wxString execute = "tetgen -" + prj.model.sld.tetgen_command + " " + prj.full_path_name + ".poly";
        std::cout << "tetgen -" << prj.model.sld.tetgen_command << " " << prj.name << ".poly" << "\n";
        wxArrayString output;
        wxExecute(execute, output, wxEXEC_SYNC);
        int size  = output.GetCount();
        for(int i=0; i < size; i++)
            std::cout << output[i] << std::endl;
        prj.task = project::RUN_TETGEN;
        DoStartThread();
    }
}

void MyFrame::setup_triangle(wxCommandEvent& WXUNUSED(event)) {
    wxTextEntryDialog* setupTriangleDialog = new wxTextEntryDialog(this,
            "Command line options:\nprq__a__uAcDjevngBPNEIOXzo_YS__iFlsCQVh", "Setup Triangle");
    setupTriangleDialog->SetValue(prj.model.sld.triangle_command);
    if (prj.model.sld.dim == 2 && setupTriangleDialog->ShowModal() == wxID_OK) {
        prj.model.sld.triangle_command = setupTriangleDialog->GetValue();
        wxString execute = "triangle -" + prj.model.sld.triangle_command + " " + prj.full_path_name + ".poly";
        std::cout << "triangle -" << prj.model.sld.triangle_command << " " << prj.name + ".poly" << "\n";
        wxArrayString output;
        wxExecute(execute, output, wxEXEC_SYNC);
        int size  = output.GetCount();
        for(int i=0; i < size; i++)
            std::cout << output[i] << std::endl;
        prj.task = project::RUN_TRIANGLE;
        DoStartThread();
    }
}

void MyFrame::refine_homogeneously(wxCommandEvent& WXUNUSED(event)) {
    prj.task = project::REFINE_HOMOGENEOUSLY;
    DoStartThread();
}

void MyFrame::run_analysis(wxCommandEvent& WXUNUSED(event)) {
    prj.task = project::ANALYZE;
    DoStartThread();
}

void MyFrame::OnImportPoly(wxCommandEvent& WXUNUSED(event)) {
    wxFileDialog* openFileDialog = new wxFileDialog(this, _("Open Triangle/Tetgen file"), "", "",
            "Triangle/Tetgen files (*.poly)|*.poly", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
    if (openFileDialog->ShowModal() == wxID_OK) {
        CurrentDocPath = openFileDialog->GetDirectory();
        wxFileName name(openFileDialog->GetFilename());
        bc_setting.clear();
        mtrl_setting.clear();
        prj.model.clear();
        prj.name = name.GetName();
        prj.data_path = openFileDialog->GetDirectory().Append(wxFILE_SEP_PATH);
        prj.full_path_name = prj.data_path + prj.name;
        m_treeCtrl->SetItemText(tree_project, std::string("Project: \"" + prj.name + "\"").data());
        prj.task = project::LOAD_POLY;
        DoStartThread();
    }
}

void MyFrame::OnImportHfss(wxCommandEvent& WXUNUSED(event)) {
    wxFileDialog* openFileDialog = new wxFileDialog(this, _("Open HFSS file"), "", "",
            "HFSS files (*.hfss)|*.hfss", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
    if (openFileDialog->ShowModal() == wxID_OK) {
        CurrentDocPath = openFileDialog->GetDirectory();
        wxFileName name(openFileDialog->GetFilename());
        bc_setting.clear();
        mtrl_setting.clear();
        prj.model.clear();
        prj.name = name.GetName();
        prj.data_path = openFileDialog->GetDirectory().Append(wxFILE_SEP_PATH);
        prj.full_path_name = prj.data_path + prj.name;
        wxDir dir(std::string(prj.full_path_name + ".hfssresults").data());
        wxString dirAddress = dir.GetName();
        wxArrayString dirList;
        dir.GetAllFiles(dirAddress, &dirList, "current.pnt", wxDIR_DIRS|wxDIR_FILES);
        if(dirList.size()) {
            if(dirList.size() > 1)
                std::cout << "More than one mesh available\n";
            wxFileName current_pnt(dirList[0]);
            prj.aux_path = current_pnt.GetPath(wxPATH_GET_VOLUME|wxPATH_GET_SEPARATOR);
        }
        m_treeCtrl->SetItemText(tree_project, std::string("Project: \"" + prj.name + "\"").data());
        prj.task = project::LOAD_HFSS;
        DoStartThread();
    }
}

void MyFrame::OnImportAedt(wxCommandEvent& WXUNUSED(event)) {
    wxFileDialog* openFileDialog = new wxFileDialog(this, _("Open AEDT file"), "", "",
            "AEDT files (*.aedt)|*.aedt", wxFD_OPEN|wxFD_FILE_MUST_EXIST);
    if (openFileDialog->ShowModal() == wxID_OK) {
        CurrentDocPath = openFileDialog->GetDirectory();
        wxFileName name(openFileDialog->GetFilename());
        bc_setting.clear();
        mtrl_setting.clear();
        prj.model.clear();
        prj.name = name.GetName();
        prj.data_path = openFileDialog->GetDirectory().Append(wxFILE_SEP_PATH);
        prj.full_path_name = prj.data_path + prj.name;
        wxDir dir(std::string(prj.full_path_name + ".aedtresults").data());
        wxString dirAddress = dir.GetName();
        wxArrayString dirList;
        dir.GetAllFiles(dirAddress, &dirList, "current.ngmesh", wxDIR_DIRS|wxDIR_FILES);
        if(dirList.size()) {
            if(dirList.size() > 1)
                std::cout << "More than one mesh available\n";
            wxFileName current_ngmesh(dirList[0]);
            prj.aux_path = current_ngmesh.GetPath(wxPATH_GET_VOLUME|wxPATH_GET_SEPARATOR);
        }
        m_treeCtrl->SetItemText(tree_project, std::string("Project: \"" + prj.name + "\"").data());
        prj.task = project::LOAD_AEDT;
        DoStartThread();
    }
}

void MyFrame::view_sld() {
    CleanVTKFrame();
    if(prj.model.sld.nodes.size() != 0) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(size_t i = 0; i< prj.model.sld.nodes.size(); i++) {
            points->InsertNextPoint(prj.model.sld.nodes[i][0], prj.model.sld.nodes[i][1], prj.model.sld.nodes[i][2]);
        }
        double geom_dim = prj.model.sld.get_geom_dim();
        pAxesActor->SetTotalLength(geom_dim,geom_dim,geom_dim);
        std::cout << "Axes dimensions = " << geom_dim << " m\n";
        // 2d
        vtkSmartPointer<vtkCellArray> edges =  vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkDoubleArray> edges_bc =  vtkSmartPointer<vtkDoubleArray>::New();
        edges_bc->SetNumberOfComponents(1);
        edges_bc->SetName("Boundaries");
        for(size_t i = 0; i < prj.model.sld.edges.size(); i++) {
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            for(size_t j = 0; j < prj.model.sld.edges[i].size(); j++) {
                line->GetPointIds()->SetId(j, prj.model.sld.edges[i][j]);
            }
            edges->InsertNextCell(line);
            edges_bc->InsertNextValue(double(prj.model.sld.edges_marker[i][0])/prj.model.sld.max_edges_marker);
        }
        // 3d
        vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkDoubleArray> faces_bc =  vtkSmartPointer<vtkDoubleArray>::New();
        faces_bc->SetNumberOfComponents(1);
        faces_bc->SetName("Boundaries");
        for(size_t i = 0; i < prj.model.sld.faces.size(); i++) {
            for(size_t j = 0; j < prj.model.sld.faces[i].polygons.size(); j++) {
                vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
                polygon->GetPointIds()->SetNumberOfIds(prj.model.sld.faces[i].polygons[j].size());
                for(size_t k = 0; k < prj.model.sld.faces[i].polygons[j].size(); k++)
                    polygon->GetPointIds()->SetId(k, prj.model.sld.faces[i].polygons[j][k]);
                faces->InsertNextCell(polygon);
                faces_bc->InsertNextValue(double(prj.model.sld.faces_marker[i][0])/prj.model.sld.max_faces_marker);
            }
        }
        vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
        polygonPolyData->SetPoints(points);
        if(prj.model.sld.dim == 2) {
            polygonPolyData->SetLines(edges);
            polygonPolyData->GetCellData()->SetScalars(edges_bc);
        } else if (prj.model.sld.dim == 3) {
            polygonPolyData->SetPolys(faces);
            polygonPolyData->GetCellData()->SetScalars(faces_bc);
        }
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polygonPolyData);
        pSldActor = vtkSmartPointer<vtkActor>::New();
        pSldActor->SetMapper(mapper);
        if(prj.model.sld.dim == 2) {
            pSldActor->GetProperty()->SetRepresentationToWireframe();
        } else if (prj.model.sld.dim == 3) {
            pSldActor->GetProperty()->SetRepresentationToSurface();
            pSldActor->GetProperty()->EdgeVisibilityOn();
            pSldActor->GetProperty()->SetOpacity(current_visibility);

        }
        pRenderer->AddActor(pSldActor);
        pRenderer->GetActiveCamera()->SetPosition(0,0,1);
        pRenderer->GetActiveCamera()->SetFocalPoint(0,0,0);
        pRenderer->GetActiveCamera()->SetViewUp(0,1,0);
        prj.model.sld.get_bounding_info();
        std::vector<std::vector<double> > bbox = prj.model.sld.bounding_box;
        pRenderer->ResetCamera(bbox[0][0], bbox[1][0],
                               bbox[0][1], bbox[1][1],
                               bbox[0][2], bbox[1][2]);
    } else if(prj.model.msh.nod_pos.size()!=0) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        for(size_t i = 0; i< prj.model.msh.nod_pos.size(); i++)
            points->InsertNextPoint(prj.model.msh.nod_pos[i][0], prj.model.msh.nod_pos[i][1], prj.model.msh.nod_pos[i][2]);
        vtkSmartPointer<vtkCellArray> poly =  vtkSmartPointer<vtkCellArray>::New();
        vtkSmartPointer<vtkDoubleArray> poly_mark =  vtkSmartPointer<vtkDoubleArray>::New();
        poly_mark->SetNumberOfComponents(1);
        poly_mark->SetName("Marker");
        for(size_t i = 0; i < prj.model.msh.fac_nodes.size(); i++) {
            if(prj.model.msh.fac_lab[i] != 0 && std::find(bcs_to_view.begin(),bcs_to_view.end(), prj.model.msh.fac_lab[i]) != bcs_to_view.end()) {
                vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
                polygon->GetPointIds()->SetNumberOfIds(3);
                for(size_t j = 0; j < 3; j++) {
                    polygon->GetPointIds()->SetId(j, prj.model.msh.fac_nodes[i][j]);
                }
                poly->InsertNextCell(polygon);
                poly_mark->InsertNextValue(double(prj.model.msh.fac_lab[i])/prj.model.msh.max_fac_marker);
            }
        }
        vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
        polygonPolyData->SetPoints(points);
        polygonPolyData->SetPolys(poly);
        polygonPolyData->GetCellData()->SetScalars(poly_mark);
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(polygonPolyData);
        pSldActor = vtkSmartPointer<vtkActor>::New();
        pSldActor->SetMapper(mapper);
        pSldActor->GetProperty()->SetRepresentationToSurface();
        pSldActor->GetProperty()->EdgeVisibilityOn();
        pSldActor->GetProperty()->SetOpacity(current_visibility);
        pRenderer->AddActor(pSldActor);
        pRenderer->GetActiveCamera()->SetPosition(0,0,1);
        pRenderer->GetActiveCamera()->SetFocalPoint(0,0,0);
        pRenderer->GetActiveCamera()->SetViewUp(0,1,0);
        pRenderer->ResetCamera();
        double geom_dim = prj.model.msh.get_geom_dim();
        pAxesActor->SetTotalLength(geom_dim,geom_dim,geom_dim);
        std::cout << "Axes dimensions = " << geom_dim << " m\n";
        prj.model.msh.get_bounding_info();
        std::vector<std::vector<double> > bbox = prj.model.msh.bounding_box;
        pRenderer->ResetCamera(bbox[0][0], bbox[1][0],
                               bbox[0][1], bbox[1][1],
                               bbox[0][2], bbox[1][2]);
    }
    m_pVTKWindow->Render();
}

void MyFrame::view_msh() {
    CleanVTKFrame();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(size_t i = 0; i< prj.model.msh.nod_pos.size(); i++)
        points->InsertNextPoint(prj.model.msh.nod_pos[i][0], prj.model.msh.nod_pos[i][1], prj.model.msh.nod_pos[i][2]);
    vtkSmartPointer<vtkCellArray> poly =  vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> poly_mark =  vtkSmartPointer<vtkDoubleArray>::New();
    poly_mark->SetNumberOfComponents(1);
    poly_mark->SetName("Marker");
    for(size_t i = 0; i < prj.model.msh.edg_nodes.size(); i++) {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        for(size_t j = 0; j < prj.model.msh.edg_nodes[i].size(); j++) {
            line->GetPointIds()->SetId(j, prj.model.msh.edg_nodes[i][j]);
        }
        poly->InsertNextCell(line);
        poly_mark->InsertNextValue(double(prj.model.msh.edg_lab[i])/prj.model.msh.max_edg_marker);
    }
    vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
    polygonPolyData->SetPoints(points);
    polygonPolyData->SetLines(poly);
    //polygonPolyData->GetCellData()->SetScalars(poly_mark);
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polygonPolyData);
    pMshActor = vtkSmartPointer<vtkActor>::New();
    pMshActor->SetMapper(mapper);
    pMshActor->GetProperty()->SetRepresentationToWireframe();
    pMshActor->GetProperty()->SetOpacity(current_visibility);
    pRenderer->AddActor(pMshActor);
    if(pSldActor != NULL) {
        pRenderer->AddActor(pSldActor);
        pSldActor->GetProperty()->SetOpacity(0.5);
        current_visibility = 1.0;
    }
    pRenderer->GetActiveCamera()->SetPosition(0,0,1);
    pRenderer->GetActiveCamera()->SetFocalPoint(0,0,0);
    pRenderer->GetActiveCamera()->SetViewUp(0,1,0);
    pRenderer->ResetCamera();
    if(prj.model.sld.nodes.size() != 0) {
        double geom_dim = prj.model.sld.get_geom_dim();
        pAxesActor->SetTotalLength(geom_dim,geom_dim,geom_dim);
        std::cout << "Axes dimensions = " << geom_dim << " m\n";
        prj.model.sld.get_bounding_info();
        std::vector<std::vector<double> > bbox = prj.model.sld.bounding_box;
        pRenderer->ResetCamera(bbox[0][0], bbox[1][0],
                               bbox[0][1], bbox[1][1],
                               bbox[0][2], bbox[1][2]);
    } else {
        double geom_dim = prj.model.msh.get_geom_dim();
        pAxesActor->SetTotalLength(geom_dim,geom_dim,geom_dim);
        std::cout << "Axes dimensions = " << geom_dim << " m\n";
        prj.model.msh.get_bounding_info();
        std::vector<std::vector<double> > bbox = prj.model.msh.bounding_box;
        pRenderer->ResetCamera(bbox[0][0], bbox[1][0],
                               bbox[0][1], bbox[1][1],
                               bbox[0][2], bbox[1][2]);
    }
    m_pVTKWindow->Render();
}

void MyFrame::CleanVTKFrame() {
    if(pSldActor != NULL)
        pRenderer->RemoveActor(pSldActor);
    if(pMshActor != NULL)
        pRenderer->RemoveActor(pMshActor);
    if(pFldActor != NULL)
        pRenderer->RemoveActor(pFldActor);
    if(scalarBar != NULL)
        pRenderer->RemoveActor(scalarBar);
    m_pVTKWindow->Render();
}

void MyFrame::DoStartThread() {
    m_pThread = new MyThread(this);
    if ( m_pThread->Run() != wxTHREAD_NO_ERROR ) {
        wxLogError("Can't create the thread!");
        delete m_pThread;
        m_pThread = NULL;
    }
}

wxThread::ExitCode MyThread::Entry() {
//    wxProgressDialog* dialog = new wxProgressDialog(wxT("Wait..."), wxT("Keep waiting..."), max, this, wxPD_AUTO_HIDE | wxPD_APP_MODAL);
    t.tic();
    m_pHandler->prj.execute_task();
    std::cout << "Elapsed " << t.strtoc() << "\n";
    wxQueueEvent(m_pHandler, new wxThreadEvent(wxEVT_COMMAND_MYTHREAD_COMPLETED));
    return (wxThread::ExitCode)0;     // success
}
MyThread::~MyThread() {
    wxCriticalSectionLocker enter(m_pHandler->m_pThreadCS);
    // the thread is being destroyed; make sure not to leave dangling pointers around
    m_pHandler->m_pThread = NULL;
}
void MyFrame::OnThreadCompletion(wxCommandEvent&) {
    if (prj.task == project::LOAD_POLY ||
            prj.task == project::LOAD_FES ||
            prj.task == project::LOAD_HFSS ||
            prj.task == project::LOAD_AEDT) {
        set_model();
        view_sld();
//        if(prj.model.msh.n_nodes > 0)
//            view_msh();
    } else if (prj.task == project::RUN_TETGEN ||
               prj.task == project::RUN_TRIANGLE ||
               prj.task == project::REFINE_HOMOGENEOUSLY) {
        std::cout << "Nodes  = " << prj.model.msh.n_nodes << "\n"
                  << "Edges  = " << prj.model.msh.n_edges << "\n"
                  << "Triangles  = " << prj.model.msh.n_faces << "\n"
                  << "Tetrahedra = " << prj.model.msh.n_tetras << "\n";
        view_msh();
    }
    prj.task = project::NONE;
    UpdateSystemMemory();
}
void MyFrame::OnThreadUpdate(wxCommandEvent&) {
    wxMessageOutputDebug().Printf("MYFRAME: MyThread update...\n");
    std::cout << "update\n";
}
void MyFrame::DoPauseThread() {
    // anytime we access the m_pThread pointer we must ensure that it won't
    // be modified in the meanwhile; since only a single thread may be
    // inside a given critical section at a given time, the following code
    // is safe:
    wxCriticalSectionLocker enter(m_pThreadCS);
    if (m_pThread) {       // does the thread still exist?
        // without a critical section, once reached this point it may happen
        // that the OS scheduler gives control to the MyThread::Entry() function,
        // which in turn may return (because it completes its work) making
        // invalid the m_pThread pointer
        if (m_pThread->Pause() != wxTHREAD_NO_ERROR )
            wxLogError("Can't pause the thread!");
    }
}
void MyFrame::OnClose(wxCloseEvent&) {
    {
        wxCriticalSectionLocker enter(m_pThreadCS);
        if (m_pThread) {
            wxMessageOutputDebug().Printf("MYFRAME: deleting thread");
            if (m_pThread->Delete() != wxTHREAD_NO_ERROR )
                wxLogError("Can't delete the thread!");
        }
    }
    while (1) {
        {
            wxCriticalSectionLocker enter(m_pThreadCS);
            if (!m_pThread) break;
        }
        wxThread::This()->Sleep(1);
    }
    Destroy();
}


#define TREE_EVENT_HANDLER(name)                                 \
void MyFrame::name(wxTreeEvent& event)                           \
{                                                                \
    wxTreeItemId item = event.GetItem();                        \
}
//               std::cout << m_treeCtrl->GetItemText(item).c_str() << " " << std::string(#name) << "\n";           \

TREE_EVENT_HANDLER(OnBeginDrag)
TREE_EVENT_HANDLER(OnBeginRDrag)
TREE_EVENT_HANDLER(OnDeleteItem)
TREE_EVENT_HANDLER(OnGetInfo)
TREE_EVENT_HANDLER(OnSetInfo)
TREE_EVENT_HANDLER(OnItemExpanded)
TREE_EVENT_HANDLER(OnItemExpanding)
TREE_EVENT_HANDLER(OnItemCollapsed)
TREE_EVENT_HANDLER(OnSelChanged)
TREE_EVENT_HANDLER(OnSelChanging)
TREE_EVENT_HANDLER(OnEndDrag)
TREE_EVENT_HANDLER(OnBeginLabelEdit)
TREE_EVENT_HANDLER(OnEndLabelEdit)
TREE_EVENT_HANDLER(OnItemCollapsing)
TREE_EVENT_HANDLER(OnTreeKeyDown)
//TREE_EVENT_HANDLER(OnItemActivated)
TREE_EVENT_HANDLER(OnItemStateClick)
TREE_EVENT_HANDLER(OnItemMenu)
//TREE_EVENT_HANDLER(OnItemRClick)

wxDEFINE_EVENT(wxEVT_COMMAND_MYTHREAD_COMPLETED, wxThreadEvent);
wxDEFINE_EVENT(wxEVT_COMMAND_MYTHREAD_UPDATE, wxThreadEvent);

wxBEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(FES_Quit, MyFrame::OnQuit)
    EVT_MENU(FES_About, MyFrame::OnAbout)
    EVT_MENU(FES_Open, MyFrame::OnOpen)
    EVT_MENU(FES_Save, MyFrame::OnSave)
    EVT_MENU(FES_Import_poly, MyFrame::OnImportPoly)
    EVT_MENU(FES_Import_hfss, MyFrame::OnImportHfss)
    EVT_MENU(FES_Import_aedt, MyFrame::OnImportAedt)
    EVT_MENU(FES_SetBGColor, MyFrame::OnSetBackgroundColor)
    EVT_MENU(FES_Run, MyFrame::run_analysis)
    EVT_MENU(FES_Tetgen, MyFrame::setup_tetgen)
    EVT_MENU(FES_Triangle, MyFrame::setup_triangle)
    EVT_MENU(FES_RefineHomogeneously, MyFrame::refine_homogeneously)
    EVT_SIZE(MyFrame::OnSize)
    EVT_SASH_DRAGGED_RANGE(MY_WINDOW_LEFT, MY_WINDOW_RIGHT, MyFrame::OnSashDrag)
    EVT_SASH_DRAGGED_RANGE(MY_WINDOW_BOTTOM, MY_WINDOW_BOTTOM, MyFrame::OnSashDrag)
    EVT_TREE_BEGIN_DRAG(FES_Tree_Ctrl, MyFrame::OnBeginDrag)
    EVT_TREE_BEGIN_RDRAG(FES_Tree_Ctrl, MyFrame::OnBeginRDrag)
    EVT_TREE_END_DRAG(FES_Tree_Ctrl, MyFrame::OnEndDrag)
    EVT_TREE_BEGIN_LABEL_EDIT(FES_Tree_Ctrl, MyFrame::OnBeginLabelEdit)
    EVT_TREE_END_LABEL_EDIT(FES_Tree_Ctrl, MyFrame::OnEndLabelEdit)
    EVT_TREE_DELETE_ITEM(FES_Tree_Ctrl, MyFrame::OnDeleteItem)
    EVT_TREE_SET_INFO(FES_Tree_Ctrl, MyFrame::OnSetInfo)
    EVT_TREE_ITEM_EXPANDED(FES_Tree_Ctrl, MyFrame::OnItemExpanded)
    EVT_TREE_ITEM_EXPANDING(FES_Tree_Ctrl, MyFrame::OnItemExpanding)
    EVT_TREE_ITEM_COLLAPSED(FES_Tree_Ctrl, MyFrame::OnItemCollapsed)
    EVT_TREE_ITEM_COLLAPSING(FES_Tree_Ctrl, MyFrame::OnItemCollapsing)
    EVT_TREE_SEL_CHANGED(FES_Tree_Ctrl, MyFrame::OnSelChanged)
    EVT_TREE_SEL_CHANGING(FES_Tree_Ctrl, MyFrame::OnSelChanging)
    EVT_TREE_KEY_DOWN(FES_Tree_Ctrl, MyFrame::OnTreeKeyDown)
    EVT_TREE_ITEM_ACTIVATED(FES_Tree_Ctrl, MyFrame::OnItemActivated)
    EVT_TREE_STATE_IMAGE_CLICK(FES_Tree_Ctrl, MyFrame::OnItemStateClick)
    EVT_TREE_ITEM_MENU(FES_Tree_Ctrl, MyFrame::OnItemMenu)
    EVT_TREE_ITEM_RIGHT_CLICK(FES_Tree_Ctrl, MyFrame::OnItemRClick)
    EVT_IDLE(MyFrame::OnIdleEvent)
    EVT_CLOSE(MyFrame::OnClose)
    EVT_COMMAND(wxID_ANY, wxEVT_COMMAND_MYTHREAD_UPDATE, MyFrame::OnThreadUpdate)
    EVT_COMMAND(wxID_ANY, wxEVT_COMMAND_MYTHREAD_COMPLETED, MyFrame::OnThreadCompletion)
wxEND_EVENT_TABLE()

IMPLEMENT_APP(MyApp)

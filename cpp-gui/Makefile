.SUFFIXES: .cpp .o

CC=g++

BINDIR= ./bin
SRCDIR= ./src
OBJDIR= ./obj

INCDIR = -I./src -I./dep/include -I./dep/include/vtk -I./dep/lib/mswu
LIBDIR = -L./dep/lib/

LIBWIN = -lkernel32 -luser32 -lgdi32 -lwinspool -lcomdlg32 -ladvapi32 -lshell32 -lole32 -loleaut32 -luuid -lcomctl32 -lwsock32 -lodbc32 -lshlwapi -lpsapi -liphlpapi -lversion -luxtheme -loleacc -lopengl32 -mwindows
LIBALG = -lsmumps -ldmumps -lcmumps -lzmumps -lmumps_common -lmpiseq -lpord -lopenblas -larpack -lgfortran -lquadmath
LIBWX = -lpthread -lwxmsw31u -liconv -lwxjpeg -lwxpng -lwxexpat -lwxregexu -lwxscintilla -lwxzlib -lwxtiff
LIBVTKALL = -lvtkChartsCore -lvtkCommonColor -lvtkCommonComputationalGeometry -lvtkCommonCore -lvtkCommonDataModel -lvtkCommonExecutionModel -lvtkCommonMath -lvtkCommonMisc -lvtkCommonSystem -lvtkCommonTransforms -lvtkDICOMParser -lvtkDomainsChemistry -lvtkdoubleconversion -lvtkexodusII -lvtkexpat -lvtkFiltersAMR -lvtkFiltersCore -lvtkFiltersExtraction -lvtkFiltersFlowPaths -lvtkFiltersGeneral -lvtkFiltersGeneric -lvtkFiltersGeometry -lvtkFiltersHybrid -lvtkFiltersHyperTree -lvtkFiltersImaging -lvtkFiltersModeling -lvtkFiltersParallel -lvtkFiltersParallelImaging -lvtkFiltersPoints -lvtkFiltersProgrammable -lvtkFiltersSelection -lvtkFiltersSMP -lvtkFiltersSources -lvtkFiltersStatistics -lvtkFiltersTexture -lvtkFiltersTopology -lvtkFiltersVerdict -lvtkfreetype -lvtkGeovisCore -lvtkgl2ps -lvtkglew -lvtkhdf5 -lvtkhdf5_hl -lvtkImagingColor -lvtkImagingCore -lvtkImagingFourier -lvtkImagingGeneral -lvtkImagingHybrid -lvtkImagingMath -lvtkImagingMorphological -lvtkImagingSources -lvtkImagingStatistics -lvtkImagingStencil -lvtkInfovisCore -lvtkInfovisLayout -lvtkInteractionImage -lvtkInteractionStyle -lvtkInteractionWidgets -lvtkIOAMR -lvtkIOAsynchronous -lvtkIOCityGML -lvtkIOCore -lvtkIOEnSight -lvtkIOExodus -lvtkIOExport -lvtkIOExportGL2PS -lvtkIOExportPDF -lvtkIOGeometry -lvtkIOImage -lvtkIOImport -lvtkIOInfovis -lvtkIOLegacy -lvtkIOLSDyna -lvtkIOMINC -lvtkIOMotionFX -lvtkIOMovie -lvtkIONetCDF -lvtkIOOggTheora -lvtkIOParallel -lvtkIOParallelXML -lvtkIOPLY -lvtkIOSegY -lvtkIOSQL -lvtkIOTecplotTable -lvtkIOVeraOut -lvtkIOVideo -lvtkIOXML -lvtkIOXMLParser -lvtkjpeg -lvtkjsoncpp -lvtklibharu -lvtklibproj -lvtklibxml2 -lvtkloguru -lvtklz4 -lvtklzma -lvtkmetaio -lvtknetcdf -lvtkogg -lvtkParallelCore -lvtkParallelDIY -lvtkpng -lvtkpugixml -lvtkRenderingAnnotation -lvtkRenderingContext2D -lvtkRenderingCore -lvtkRenderingFreeType -lvtkRenderingGL2PSOpenGL2 -lvtkRenderingImage -lvtkRenderingLabel -lvtkRenderingLOD -lvtkRenderingOpenGL2 -lvtkRenderingSceneGraph -lvtkRenderingUI -lvtkRenderingVolume -lvtkRenderingVolumeOpenGL2 -lvtkRenderingVtkJS -lvtksqlite -lvtksys -lvtkTestingRendering -lvtktheora -lvtktiff -lvtkverdict -lvtkViewsContext2D -lvtkViewsCore -lvtkViewsInfovis -lvtkWrappingTools -lvtkzlib
LIBVTK = -lvtkRenderingOpenGL2 -lvtkRenderingVolumeOpenGL2 -lvtkRenderingVolume -lvtkRenderingUI -lvtkRenderingFreeType -lvtkRenderingAnnotation -lvtkInteractionStyle -lvtkRenderingCore -lvtkCommonColor -lvtkFiltersCore -lvtkFiltersSources -lvtkFiltersGeneral -lvtkImagingCore -lvtkCommonExecutionModel -lvtkCommonDataModel -lvtkCommonMath -lvtkCommonTransforms -lvtkCommonSystem -lvtkCommonCore -lvtkCommonMisc -lvtksys -lvtkloguru -lvtkglew -lvtkfreetype -lvtkzlib

CFLAGS = $(INCDIR) -std=gnu++11 -m64 -O2 -static -fopenmp -DNDEBUG -D_FILE_OFFSET_BITS=64 -DWX_PRECOMP -fno-common -fpermissive -Wall -Wundef -Wunused-parameter -Wno-ctor-dtor-privacy -Woverloaded-virtual -Wno-deprecated-declarations
LFLAGS = $(LIBDIR) -m64 -fopenmp -static -s $(LIBALG) $(LIBWX) $(LIBVTK) $(LIBWIN)

OBJS = $(addprefix $(OBJDIR)/, wxfes.o model.o project.o solver.o wxVTKRenderWindowInteractor.o)

ifdef OS
 RM = del /F /S /Q
 FixPath = $(subst /,\,$1)
else
 ifeq ($(shell uname), Linux)
  RM = rm -f
  FixPath = $1
 endif
endif

all: $(OBJS)
	$(CC) -o $(BINDIR)/wxfes $(OBJS) $(LFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(call FixPath,$(OBJDIR)/*.o)



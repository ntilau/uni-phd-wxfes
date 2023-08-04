.SUFFIXES: .cpp .o

CC=gcc
CCP=g++

BINDIR= ..
SRCDIR= .
OBJDIR= ./obj

INCDIR = -I.
LIBDIR-WIN = -L.
LIBDIR-LNX = -L.

ifdef OS
   RM = del /F /S /Q
   FixPath = $(subst /,\,$1)
   LIBDIR = $(LIBDIR-WIN)
else
   ifeq ($(shell uname), Linux)
      RM = rm -f
      FixPath = $1
      LIBDIR = $(LIBDIR-LNX)
   endif
endif

CFLAGS = $(INCDIR) -Wall -m64 -O2
CPPFLAGS = $(INCDIR) -m64 -O2 -fpermissive
LFLAGS = $(LIBDIR) -m64 -s -lm
 
OBJS = $(addprefix $(OBJDIR)/, triangle.o matfiles.o main.o)

all: $(OBJS)
	$(CC) -o $(BINDIR)/triangle.exe $(addprefix $(OBJDIR)/, triangle.o) $(LFLAGS)
	$(CC) -o $(BINDIR)/IOrMesh.exe $(addprefix $(OBJDIR)/, matfiles.o main.o ) $(LFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c  $< -o $@

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CCP) $(CPPFLAGS) -c $< -o $@

.PHONY: clean
clean:
	$(RM) $(call FixPath,$(OBJDIR)/*.o)

.RECIPEPREFIX := >
PROGRAM := rotex
BINDIR ?= bin
PRGDIR ?= main
SRCDIR ?= src
OBJDIR ?= build/obj
MODDIR ?= build/mod
SRCDIR_FACE ?= external/FACE/src/lib
SRCDIR_FORBEAR ?= external/forbear/src/lib
SRCDIR_BSPLINE_FORTRAN  ?= external/bspline-fortran/src
SRCDIR_CDMSREADER ?= external/CDMSreader/src
SRCDIR_WIGNERD ?= external/WignerD/src
VPATH := $(SRCDIR):$(SRCDIR_FACE):$(SRCDIR_FORBEAR):$(SRCDIR_CDMSREADER):$(SRCDIR_BSPLINE_FORTRAN):$(SRCDIR_WIGNERD):$(PRGDIR)
USE_FORBEAR    ?= 1
USE_OPENMP     ?= 1
USE_CDMSREADER ?= 1
FC := gfortran
FFLAGS ?= -cpp \
          -I $(MODDIR) \
          -J $(MODDIR) \
          -D_REENTRANT \
          -mtune=generic \
          -march=native \
          -O3 \
					-g -fbacktrace \
          -Wimplicit-interface \
          -Werror=implicit-interface \
          -fmax-errors=1 \
          -fcoarray=single \
          -funroll-loops \
          -fPIC \
          -fimplicit-none \
          -ffree-form
FLFLAGS ?= -l lapack

ifeq ($(USE_OPENMP),1)
  FFLAGS += -fopenmp -DUSE_OPENMP
endif

# -- dependencies
ifeq ($(USE_FORBEAR),1)
  FFLAGS += -DUSE_FORBEAR
  SRC_FACE := $(shell find $(SRCDIR_FACE) -type f \( -name '*.f' -o -name '*.f90' -o -name '*.F90' \))
  OBJ_FACE := $(notdir $(SRC_FACE))
  OBJ_FACE := $(OBJ_FACE:.f90=.o)
  OBJ_FACE := $(OBJ_FACE:.F90=.o)
  OBJ_FACE := $(OBJ_FACE:.f=.o)
  OBJ_FACE := $(addprefix $(OBJDIR)/,$(OBJ_FACE))

  SRC_FORBEAR := $(shell find $(SRCDIR_FORBEAR) -type f \( -name '*.f' -o -name '*.f90' -o -name '*.F90' \))
  OBJ_FORBEAR := $(notdir $(SRC_FORBEAR))
  OBJ_FORBEAR := $(OBJ_FORBEAR:.f90=.o)
  OBJ_FORBEAR := $(OBJ_FORBEAR:.F90=.o)
  OBJ_FORBEAR := $(OBJ_FORBEAR:.f=.o)
  OBJ_FORBEAR := $(addprefix $(OBJDIR)/,$(OBJ_FORBEAR))
endif
ifeq ($(USE_CDMSREADER),1)
  FFLAGS += -DUSE_CDMSREADER
  SRC_CDMSREADER := $(shell find $(SRCDIR_CDMSREADER) -type f \( -name '*.f' -o -name '*.f90' -o -name '*.F90' \))
  OBJ_CDMSREADER := $(notdir $(SRC_CDMSREADER))
  OBJ_CDMSREADER := $(OBJ_CDMSREADER:.f90=.o)
  OBJ_CDMSREADER := $(OBJ_CDMSREADER:.F90=.o)
  OBJ_CDMSREADER := $(OBJ_CDMSREADER:.f=.o)
  OBJ_CDMSREADER := $(addprefix $(OBJDIR)/,$(OBJ_CDMSREADER))
endif
# -- Bspline Fortran
SRC_BSPLINE_FORTRAN := $(shell find $(SRCDIR_BSPLINE_FORTRAN) -type f \( -name '*.f' -o -name '*.f90' -o -name '*.F90' \))
OBJ_BSPLINE_FORTRAN := $(notdir $(SRC_BSPLINE_FORTRAN))
OBJ_BSPLINE_FORTRAN := $(OBJ_BSPLINE_FORTRAN:.f90=.o)
OBJ_BSPLINE_FORTRAN := $(OBJ_BSPLINE_FORTRAN:.F90=.o)
OBJ_BSPLINE_FORTRAN := $(OBJ_BSPLINE_FORTRAN:.f=.o)
OBJ_BSPLINE_FORTRAN := $(addprefix $(OBJDIR)/,$(OBJ_BSPLINE_FORTRAN))
# -- WignerD
SRC_WIGNERD:= $(shell find $(SRCDIR_WIGNERD) -type f \( -name '*.f' \))
OBJ_WIGNERD := $(notdir $(SRC_WIGNERD))
OBJ_WIGNERD := $(OBJ_WIGNERD:.f=.o)
OBJ_WIGNERD := $(addprefix $(OBJDIR)/,$(OBJ_WIGNERD))

# -- ROTEX src and object files
PRG := $(shell find $(PRGDIR) -name '*.f')
SRC := $(shell find $(SRCDIR) -name '*.f')
OBJ := $(notdir $(SRC))
OBJ := $(OBJ:.f=.o)
OBJ := $(addprefix $(OBJDIR)/,$(OBJ))

# -- all objects
OBJ_ALL := $(OBJ_FACE) $(OBJ_FORBEAR) $(OBJ_BSPLINE_FORTRAN) $(OBJ_CDMSREADER) $(OBJ_WIGNERD) $(OBJ)

# -- default target
build: dirs $(BINDIR)/$(PROGRAM)

# -- ensure existence of necessary folders
dirs:
> mkdir -p $(OBJDIR) $(MODDIR) $(BINDIR)

# -- dependency build dependencies
ifeq ($(USE_FORBEAR),1)
  $(OBJ_FORBEAR): $(OBJ_FACE)
  $(OBJDIR)/forbear_element_object.o: $(OBJDIR)/forbear_kinds.o
  $(OBJDIR)/forbear_bar_object.o: $(OBJDIR)/forbear_element_object.o
	$(OBJDIR)/forbear.o: $(OBJDIR)/forbear_bar_object.o
endif
ifeq ($(USE_CDMSREADER),1)
	$(OBJDIR)/CDMSreader__types.o: $(OBJDIR)/CDMSreader__system.o
  $(OBJDIR)/CDMSreader__constants.o: $(OBJDIR)/CDMSreader__types.o
	$(OBJDIR)/CDMSreader__readwrite.o: $(OBJDIR)/CDMSreader__constants.o $(OBJDIR)/CDMSreader__types.o $(OBJDIR)/CDMSreader__system.o
endif
# -- bspline deps
$(OBJDIR)/bspline_defc_module.o: $(OBJDIR)/bspline_kinds_module.o $(OBJDIR)/bspline_blas_module.o
$(OBJDIR)/bspline_sub_module.o: $(OBJDIR)/bspline_kinds_module.o
$(OBJDIR)/bspline_oo_module.o: $(OBJDIR)/bspline_kinds_module.o $(OBJDIR)/bspline_sub_module.o
$(OBJDIR)/bspline_module.o: $(OBJDIR)/bspline_oo_module.o $(OBJDIR)/bspline_sub_module.o $(OBJDIR)/bspline_kinds_module.o $(OBJDIR)/bspline_defc_module.o $(OBJDIR)/bspline_blas_module.o
# -- WignerD deps
$(OBJDIR)/wignerd.o: $(OBJDIR)/wignerd__types.o $(OBJDIR)/wignerd__constants.o $(OBJDIR)/wignerd__characters.o $(OBJDIR)/wignerd__system.o $(OBJDIR)/wignerd__functions.o
$(OBJDIR)/wignerd__characters.o: $(OBJDIR)/wignerd__types.o $(OBJDIR)/wignerd__constants.o
$(OBJDIR)/wignerd__functions.o: $(OBJDIR)/wignerd__types.o $(OBJDIR)/wignerd__constants.o
$(OBJDIR)/wignerd__system.o: $(OBJDIR)/wignerd__types.o
$(OBJDIR)/wignerd__constants.o: $(OBJDIR)/wignerd__types.o

# -- ROTEX build dependencies
$(OBJDIR)/rotex__arrays.o: $(OBJDIR)/rotex__types.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o
$(OBJDIR)/rotex__characters.o: $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__constants.o
$(OBJDIR)/rotex__cbxs.o: $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__types.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__wigner.o $(OBJDIR)/rotex__hypergeometric.o
$(OBJDIR)/rotex__constants.o: $(OBJDIR)/rotex__kinds.o
$(OBJDIR)/rotex__drivers.o: $(OBJDIR)/rotex__arrays.o $(OBJDIR)/rotex__cbxs.o $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__hamilton.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__mqdtxs.o $(OBJDIR)/rotex__reading.o $(OBJDIR)/rotex__rft.o $(OBJDIR)/rotex__splines.o $(OBJDIR)/rotex__symmetry.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o $(OBJDIR)/rotex__utils.o $(OBJDIR)/rotex__writing.o $(OBJ_CDMSREADER)
$(OBJDIR)/rotex__functions.o: $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__system.o
$(OBJDIR)/rotex__hamilton.o: $(OBJDIR)/rotex__arrays.o $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__linalg.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o $(OBJDIR)/rotex__utils.o
$(OBJDIR)/rotex__hypergeometric.o: $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__polygamma.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__utils.o
$(OBJDIR)/rotex__linalg.o: $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__system.o
$(OBJDIR)/rotex__mqdtxs.o: $(OBJDIR)/rotex__arrays.o $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__linalg.o $(OBJDIR)/rotex__progress.o $(OBJDIR)/rotex__symmetry.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o
$(OBJDIR)/rotex__polygamma.o: $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__utils.o
$(OBJDIR)/rotex__progress.o: $(OBJ_FORBEAR)
$(OBJDIR)/rotex__reading.o: $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__symmetry.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o $(OBJDIR)/rotex__utils.o
$(OBJDIR)/rotex__rft.o: $(OBJDIR)/rotex__arrays.o $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__linalg.o $(OBJDIR)/rotex__symmetry.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o $(OBJDIR)/rotex__wigner.o
$(OBJDIR)/rotex__splines.o: $(OBJDIR)/bspline_module.o $(OBJDIR)/rotex__arrays.o $(OBJDIR)/rotex__kinds.o
$(OBJDIR)/rotex__symmetry.o: $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o
$(OBJDIR)/rotex__types.o: $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__system.o
$(OBJDIR)/rotex__utils.o: $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o
$(OBJDIR)/rotex__wigner.o: $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o
$(OBJDIR)/rotex__writing.o: $(OBJDIR)/rotex__arrays.o $(OBJDIR)/rotex__characters.o $(OBJDIR)/rotex__constants.o $(OBJDIR)/rotex__functions.o $(OBJDIR)/rotex__kinds.o $(OBJDIR)/rotex__symmetry.o $(OBJDIR)/rotex__system.o $(OBJDIR)/rotex__types.o

# -- dependencies
$(OBJ): $(OBJ_FACE) $(OBJ_FORBEAR) $(OBJ_CDMSREADER) $(OBJ_BSPLINE_FORTRAN) $(OBJ_WIGNERD)

# -- link
$(BINDIR)/$(PROGRAM): $(OBJ_ALL) $(PRGDIR)/rotex.f | dirs
> $(FC) $(FFLAGS) $(FLFLAGS) -o $@ $(OBJ_ALL) $(PRGDIR)/rotex.f

# -- general build command for .f90, .F90, and .f files
$(OBJDIR)/%.o: %.f90 | dirs
> $(FC) $(FFLAGS) -c $< -o $@
$(OBJDIR)/%.o: %.F90 | dirs
> $(FC) $(FFLAGS) -c $< -o $@
$(OBJDIR)/%.o: %.f | dirs
> $(FC) $(FFLAGS) -c $< -o $@

clean:
>	rm -f $(MODDIR)/* $(OBJDIR)/* $(BINDIR)/*

# -- for when stuff doesn't work
debug:
> @echo ""
> @echo "---------- debug ----------"
> @echo ""
> @echo " FC = $(FC)"
> @echo " BINDIR = $(BINDIR)"
> @echo " SRCDIR = $(SRCDIR)"
> @echo " OBJDIR = $(OBJDIR)"
> @echo " PRGDIR = $(PRGDIR)"
> @echo " MODDIR = $(MODDIR)"
> @echo " SRCDIR_FACE = $(SRCDIR_FACE)"
> @echo " SRCDIR_FORBEAR = $(SRCDIR_FORBEAR)"
> @echo " SRCDIR_BSPLINE_FORTRAN = $(SRCDIR_BSPLINE_FORTRAN) "
> @echo " SRCDIR_CDMSREADER = $(SRCDIR_CDMSREADER)"
> @echo ""
> @echo "  ROTEX source files:"
> @echo "    $(SRC)"
> @echo ""
> @echo "  USE_OPENMP: $(USE_OPENMP)"
> @echo "  USE_FORBEAR: $(USE_FORBEAR)"
> @echo "  USE_CDMSREADER: $(USE_CDMSREADER)"
> @echo ""
> @echo "  Bspline-fortran source files:"
> @echo "    $(SRC_BSPLINE_FORTRAN)"
> @echo ""
> @echo "  WignerD source files:"
> @echo "    $(SRC_WIGNERD)"
> @echo ""
> @echo "  FACE source files:"
> @echo "    $(SRC_FACE)"
> @echo ""
> @echo "  Forbear source files:"
> @echo "    $(SRC_FORBEAR)"
> @echo ""
> @echo "  CDMSreader source files:"
> @echo "    $(SRC_CDMSREADER)"
> @echo ""
> @echo "  object files:"
> @echo "    $(OBJ_ALL)"

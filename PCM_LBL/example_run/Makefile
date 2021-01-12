# shell = /bin/sh
raydir  = $(shell pwd)
un      = $(shell whoami)
host    = $(shell echo `hostname`)
exec    = PCM_LBL.e
for     = gfortran
make    = make
osv     = linux
ostypev = $(shell echo ${ostype})
warning =
debug   = -g -cpp -fbounds-check -ffpe-trap='invalid','zero' -Wall -W -Wsurprising -pedantic -O -fstrict-aliasing -g3 -ggdb3 -fno-backtrace
fflags  = $(warning) #$(debug)
inc=-I./include

# dir is the path to source files directory
# obj is the path to object files directory
# src_prog is the list of source files to be compiled
# includes is -i+path to directory that contains 'include' files
# Note: need to improve or replace ioipsl_getincom

dir = /Users/robin/Downloads/mars_redox_2020-master/PCM_LBL/source
obj = /Users/robin/Downloads/mars_redox_2020-master/PCM_LBL/objects

src_prog = \
	$(dir)/params_consts.f90 \
	$(dir)/composition.f90 \
	$(dir)/thermodynamics.f90 \
	$(dir)/atmos_profile.f90 \
	$(dir)/atmos_structure.f90 \
	$(dir)/cia.f90 \
	$(dir)/line_profs.f90 \
	$(dir)/cross_section.f90 \
	$(dir)/planck.f90 \
	$(dir)/temperature_k_grid.f90 \
	$(dir)/longwave.f90 \
	$(dir)/shortwave.f90 \
	$(dir)/radiation.f90 \
	$(dir)/main.f90 \

# what to compile
obj_prog=$(src_prog:$(dir)/%.for=$(obj)/%.o)
# how to compile - not necessary??
$(obj)/%.o: $(dir)/%.for $(for) $(fflags) $(inc) -c -o $@ $(@:$(obj)/%.o=$(dir)/%.for)

clean :
	@rm -f *.o $(obj)/*.o *.oub *.oup *.oug *.mod
	@rm -f *%
	@rm -f *~
	@rm -f $(dir)/*.~
	@rm -f $(dir)/*.*~
	@rm -f *.eps
	@rm -f *.ps
	@rm -f *.jpg
	@echo 'files have been removed'

where :
	@echo current directory: $(raydir)
	@echo source directory: $(dir)

info :
	@echo host name: $(hostname)
	@echo operating system: $(osv)
	@echo type of system: $(ostypev)
	@echo user name: $(un)

all :
	@echo $(dir)
	@echo '------ initialization -------'
	make -f Makefile PCM_LBL

PCM_LBL :	$(obj_prog)
	@echo '------ compilation -------'
	$(for) $(fflags) -o $(exec) $(obj_prog)
	@echo '------ cleanup -------'
	@rm -f *.mod
	@echo 'all .mod files removed.'
	@echo '--- use program ' $(exec)  '(compiled for' $(osv) $(ostypev) 'by' $(un) ')'

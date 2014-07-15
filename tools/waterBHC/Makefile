# Copy updated .mat files to the ../../data/ directory

orgFiles = muEff140Sn.mat muEff80.mat polycr140Sn.mat polycr80.mat
dstFiles = $(addprefix ../../data/,$(orgFiles))

all : 
	@echo Usage: make copy_mat_files

copy_mat_files : $(dstFiles)

../../data/%.mat : %.mat
	cp $< ../../data

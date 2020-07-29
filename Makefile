target = main.exe
objs=read_input.o\
     cal_grid.o\
     glo_data.o\
     alloc_field.o\
     dealloc_field.o\
     ini_J.o\
     FCT_main.o\
# linking the target depends on the objects
#typ = mpif90
typ = f90
opt = -free
$(target):$(objs)
	$(typ) $(opt) $(objs) -o $(target)
#dependencies
FCT_main.o:FCT_main.f90 read_input.o glo_data.o alloc_field.o dealloc_field.o ini_J.o
	$(typ) $(opt) -c FCT_main.f90
glo_data.o:glo_data.f90
	$(typ) $(opt) -c glo_data.f90
read_input.o:read_input.f90 glo_data.o
	$(typ) $(opt) -c read_input.f90
cal_grid.o:cal_grid.f90 glo_data.o
	$(typ) $(opt) -c cal_grid.f90
alloc_field.o:alloc_field.f90
	$(typ) $(opt) -c alloc_field.f90
dealloc_field.o:dealloc_field.f90
	$(typ) $(opt) -c dealloc_field.f90
ini_J.o:ini_J.f90 glo_data.o
	$(typ) $(opt) -c ini_J.f90
new:clean $(target)
clean:
	rm -fr $(objs)

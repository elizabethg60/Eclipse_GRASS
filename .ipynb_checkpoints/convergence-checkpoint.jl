using Revise
using MyProject

MyProject.get_kernels()

# lat_array = range(115, 115, step = 5)

# for i in lat_array
#     MyProject.gottingen_loop(i)
# end

sub_array = range(14, 16, step = 2)

for i in sub_array
    MyProject.gottingen_loop(197, i)
end
using Test, MyProject

#tests ran for max eclipse of gottingen 2015 eclipse 

epoch = 4.80116587185603e8 #utc2et("2015-03-20T09:42:00")

@testset "Test Set 1: Unit Tests" begin
	#confirm recovered weighted velocity decreases as grid size increases in the case of no moon
	@test MyProject.compute_rv(100,100, epoch, 0, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1] > MyProject.compute_rv(50,50, epoch, 0, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
	#above test for parallel code
	@test MyProject.compute_rv_pa(100,100, epoch, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1] > MyProject.compute_rv_pa(50,50, epoch, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
end

@testset "Test Set 2: Integration Tests" begin
	#confirm that if moon radius >> sun radius that recovered velocity is nan (nothing visible)
    @test isnan(MyProject.compute_rv(100,100, epoch, 0, 9.905548, 51.54548, 0.15, "optical", moon_r = 696340*5.0)[1])
	#above test for parallel code
	@test isnan(MyProject.compute_rv_pa(100,100, epoch, 9.905548, 51.54548, 0.15, "optical", moon_r = 696340*5.0)[1])

	#confirm that is no moon present then recovered velocity approx the same at max epoch and first timestamp 
	max = MyProject.compute_rv(100,100, epoch, 0, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
	first = MyProject.compute_rv(100,100, 4.8010716718560225e8, 0, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
	@test isapprox(max,first; rtol = 1)
	#above test for parallel code
	max = MyProject.compute_rv_pa(100,100, epoch, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
	first = MyProject.compute_rv_pa(100,100, 4.8010716718560225e8, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
	@test isapprox(max,first; rtol = 1)

	#integration test between serial and parallel code
	max_serial = MyProject.compute_rv(100,100, epoch, 0, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]

	max_parallel = MyProject.compute_rv_pa(100,100, epoch, 9.905548, 51.54548, 0.15, "optical", moon_r = 0.0)[1]
	@test isapprox(max_serial, max_parallel; rtol = 1)
end  
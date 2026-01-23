using TestItems
using TestItemRunner
@run_package_tests

@testitem "_" begin
    import Aqua
    Aqua.test_all(Synchray;
		undefined_exports=(;broken=true),  # reexport StaticArrays and Accessors that both define insert()
	)

    import CompatHelperLocal as CHL
    CHL.@check()
end

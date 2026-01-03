using Preferences
if get(ENV, "GITHUB_ACTIONS", "false") == "true"
	@info "Running in GitHub Actions; set dispatch_doctor_mode to error"
	set_preferences!("Synchray", "dispatch_doctor_mode" => "error")
end

using TestItems
using TestItemRunner
@run_package_tests

# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

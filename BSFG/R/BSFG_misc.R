reinstall_BSFG = function(local_library = '~/R/x86_64-pc-linux-gnu-library/3.3/')
{
  BSFG_path = 'https://github.com/deruncie/SparseFactorMixedModel'
  withr::with_libpaths(local_library,devtools::install_git(BSFG_path,branch = 'develop',subdir = 'BSFG'))
}

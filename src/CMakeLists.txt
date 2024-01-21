add_custom_target(bio_condensates)

add_amdis_executable(NAME cahn_hilliard SOURCES cahn_hilliard.cc DIM 2 DOW 2 ALBERTA_GRID)
add_amdis_executable(NAME ch_curve SOURCES ch_curve.cc DIM 2 DOW 2 ALBERTA_GRID)
add_amdis_executable(NAME ch_active_b_plus SOURCES ch_active_b_plus.cc DIM 2 DOW 2 ALBERTA_GRID)

add_dependencies(bio_condensates
    cahn_hilliard
    ch_curve
    ch_active_b_plus
)

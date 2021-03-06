
fap_dic = {
        "just_repeat":
                   ["dft_code = vasp",
                    "dispersion = False",
                    "charge_method = repeat",
                    "no_gcmc = True",
                    "no_properties = True",
                    "optim_cell = False",
                    "optim_h = False",
                    "optim_all = False",
                    "quiet = True",
                    ],

        "vasp_opt":
                   ["dft_code = vasp",
                    "dispersion = True", 
                    "charge_method = repeat",
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False",
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "no_properties = False",
                    "optim_cell = True", 
                    "optim_all = True", 
                    "quiet = True"
                    ],

        "uff_qeq":
                   ["charge_method = gulp",
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False",
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False",
                    "optim_h = False",
                    "quiet = True",
                    "qeq_parameters =",
                    "    Cu    4.20000000    4.22000000"
                    ],

        "vasp_gcmc_sp":
                   ["charge_method = repeat",
                    "dft_code = vasp",
                    "dispersion = False",
                    "vasp_ncpu = 24",
                    "guests = CO2",
                    "mc_probability_plot = True",
                    "fold = True",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "quiet = True"
                    ],

        "noq_gcmc":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False", 
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "no_charges = True",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "quiet = True"
                    ],

        "egulp_3":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False", 
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "charge_method = egulp",
                    "qeq_parameters =", 
                    "    Cu    5.28660000    3.28889297",
                    "    Zn    3.47871363    4.14128701",
                    "     C    5.448105      5.53192500",
                    "     N    6.30432       6.03174000",
                    "     O   10.70614       8.57850000",
                    "     F    6.364800     11.136000",
                    "    Cl    5.824000      7.272000",
                    "    Br    5.760000      8.760000",
                    "     I    5.426667      5.720000",
                    "     S    3.080533      4.336800",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "quiet = True"
                    ],

        "egulp_4d":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False", 
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "charge_method = egulp",
                    "qeq_parameters =", 
                    "    Cu    5.42900000    3.46800000",
                    "    Zn    3.70100000    4.46300000",
                    "     C    5.431000      5.85700000",
                    "     N    6.68800       6.62200000",
                    "     O    8.71400       8.56800000",
                    "     F    6.416000     11.131000",
                    "    Cl    5.821000      7.273000",
                    "    Br    5.692000      8.760000",
                    "     I    5.431000      5.720000",
                    "     S    3.369000      5.092000",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "quiet = True"
                    ],


        "egulp_4dt2":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False", 
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "charge_method = egulp",
                    "qeq_parameters =", 
                    "    Cu    5.429000      3.468000",
                    "    Zn    3.701000      4.463000",
                    "     C    5.431000      5.857000",
                    "     N    6.688000      6.622000",
                    "     O    8.714000      8.568000",
                    "     F    6.416000     11.131000",
                    "    Cl    5.821000      7.273000",
                    "    Br    5.692000      8.760000",
                    "     I    5.431000      5.720000",
                    "     S    3.369000      5.092000",
                    "    800   8.714000      8.568000",
                    "    801  10.597000      9.744000",
                    "    802   7.968000     10.323000",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "egulp_typed_atoms = True",
                    "quiet = True"
                    ],

        "egulp_4dt3":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = True",
                    "fold = True",
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "charge_method = egulp",
                    "qeq_parameters =", 
                    "    Cu    5.429000      3.468000",
                    "    Zn    3.701000      4.463000",
                    "     C    5.431000      5.857000",
                    "     N    6.688000      6.622000",
                    "     O    8.714000      8.568000",
                    "     F    6.416000     11.131000",
                    "    Cl    5.821000      7.273000",
                    "    Br    5.692000      8.760000",
                    "     I    5.431000      5.720000",
                    "     S    3.369000      5.092000",
                    "    800   8.714000      8.568000",
                    "    801  10.528000      9.543000",
                    "    802   8.086000     10.187000",
                    "   1001   4.035000      6.722000",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "egulp_typed_atoms = True",
                    "quiet = True"
                    ],
        "egulp_4dt3v1":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = False", 
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "charge_method = egulp",
                    "qeq_parameters =", 
                    "    Cu    5.429000      3.468000",
                    "    Zn    3.701000      4.463000",
                    "     V    3.846000      3.667000",
                    "     C    5.431000      5.857000",
                    "     N    6.688000      6.622000",
                    "     O    8.714000      8.568000",
                    "     F    6.416000     11.131000",
                    "    Cl    5.821000      7.273000",
                    "    Br    5.692000      8.760000",
                    "     I    5.431000      5.720000",
                    "     S    3.369000      5.092000",
                    "    800   8.714000      8.568000",
                    "    801  10.528000      9.543000",
                    "    802   8.086000     10.187000",
                    "   1001   4.035000      6.722000",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "egulp_typed_atoms = True",
                    "quiet = True"
                    ],

        "egulp_4dt3v2":
                   ["find_maxima = False", 
                    "guests = CO2",
                    "mc_eq_steps = 2000000",
                    "mc_prod_steps = 10000000",
                    "mc_numguests_freq = 50000",
                    "mc_probability_plot = True",
                    "fold = True",
                    "mc_temperature = 298.0",
                    "mc_pressure = 0.15",
                    "charge_method = egulp",
                    "qeq_parameters =", 
                    "    Cu    5.429000      3.468000",
                    "    Zn    3.701000      4.463000",
                    "     V    4.093000      4.217000",
                    "     C    5.431000      5.857000",
                    "     N    6.688000      6.622000",
                    "     O    8.714000      8.568000",
                    "     F    6.416000     11.131000",
                    "    Cl    5.821000      7.273000",
                    "    Br    5.692000      8.760000",
                    "     I    5.431000      5.720000",
                    "     S    3.369000      5.092000",
                    "    800   8.714000      8.568000",
                    "    801  10.528000      9.543000",
                    "    802   8.086000     10.187000",
                    "   1001   4.035000      6.722000",
                    "no_dft = True",
                    "optim_all = False",
                    "optim_cell = False", 
                    "optim_h = False",
                    "egulp_typed_atoms = True",
                    "quiet = True"
                    ]
}

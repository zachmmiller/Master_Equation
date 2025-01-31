-- require("io") -- uncomment to print stuff using "io.write()"
dofile("/Users/zachmiller/Programming/C++/Master-Equation/lua_libraries/csv.lua")

e_min = 0 -- wavenumber
e_max = 20000
e_step = 20

t_min = 0 -- seconds
t_max = 2
t_step = 0.001

temperature = 120 + 273.15 -- kelvin

modes = csv.load(config_directory .. "/vibrations.csv", ",", true)

vib_modes = csv_col(modes, 1) -- wavenumber
ir_intens = csv_col(modes, 2) -- km/mol
vib_degen = csv_col(modes, 3) -- count

TS_vib_modes = csv_col(modes, 1)
TS_vib_degen = csv_col(modes, 3)

e0 = 7345 -- wavenumber

initially_boltzmann = true
initial_temperature = temperature
initial_integral_abundance = 1 -- Sum of all abundance

no_RRKM = false

save_modes = true
save_initial_condition = true
save_boltzmann = true
save_eigenvalues = true
save_time_data = true

output_directory = config_directory .. "/../outputs"

number_of_threads = 8

-- Scripting

TS_state_mode_index = 7
TS_low_mode_scalar = 1.05
TS_vib_modes[1] = TS_vib_modes[1] * TS_low_mode_scalar
TS_vib_modes[2] = TS_vib_modes[2] * TS_low_mode_scalar
TS_vib_degen[TS_state_mode_index] = TS_vib_degen[TS_state_mode_index]-1

vib_mode_scalar = 0.98
ir_intens_scalar = 1.0

for i, v in ipairs(vib_modes)
  do
  vib_modes[i] = v * vib_mode_scalar
  end

for i, v in ipairs(TS_vib_modes)
  do TS_vib_modes[i] = v * vib_mode_scalar
  end

for i, I in ipairs(ir_intens)
  do ir_intens[i] = I * ir_intens_scalar
  end

-- init_pop = slice_csv(csv.load(config_directory .. "/initial_population.csv", ",", false), 1)

function initial_population(index)
  energy = e_step * index + e_min
  return energy
end



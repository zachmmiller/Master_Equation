-- require("io") -- uncomment to print stuff using "io.write()"
dofile(config_directory .. "/../../../../lua_libraries/csv.lua") -- string are concatenated in Lua with the ".." operator.

energies = {10000, 25000, 50000, 100000, 150000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 2000000, 3000000, 4000000} -- wavenumber

modes = csv.load(config_directory .. "/../vibrations.csv", ",", true)

vib_modes = csv_col(modes, 1) -- wavenumber
ir_intens = csv_col(modes, 2) -- km/mol
vib_degen = csv_col(modes, 3) -- count

output_directory = config_directory .. "/../outputs"

number_of_threads = 8

-- Scripting

vib_mode_scalar = 0.98

for i, v in ipairs(vib_modes)
  do
  vib_modes[i] = v * vib_mode_scalar
  end






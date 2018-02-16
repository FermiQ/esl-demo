--[[
Example Lua script to exchange data with ESL-DEMO
--]]

function esl_comm()

   esl.receive({'Geometry.xyz'})
   ret_tbl = {}

   if esl.state == esl.INITIALIZE then
      -- Uncomment to see available variables to exchange
      -- esl.print_allowed()
   end

   if esl.state == esl.SCF_LOOP then
      ret_tbl = scf(esl)
   end

end



res = {}
local math = require "math"
function scf(esl)
   -- Start by
   esl.IOprint('# LUA: Currently in SCF!')
   --esl.print_allowed()

   -- Retrieve the residue
   esl.receive({'SCF.Residue', 'SCF.Mixer.Alpha'})

   -- Now we can access
   res[#res + 1] = esl.SCF.Residue
   print('# Lua residue: ', res[#res])

   -- Check if we can update the mixing parameter
   if #res >= 2 then
      print('# LUA: Updating mixing weight (from): ', esl.SCF.Mixer.Alpha)
      if res[#res] - res[#res-1] < 0 then
	 -- Increase by 5%
	 esl.SCF.Mixer.Alpha = esl.SCF.Mixer.Alpha * 1.05
      else
	 -- Decrease by 10%
	 esl.SCF.Mixer.Alpha = esl.SCF.Mixer.Alpha * 0.9
      end
      esl.SCF.Mixer.Alpha = math.min(0.9, esl.SCF.Mixer.Alpha)

      esl.send({'SCF.Mixer.Alpha'})
   end
end
   

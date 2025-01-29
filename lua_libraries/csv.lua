local function parse_row(input, sep, pos)
	local row = {}
	local pos = pos or 1
	--io.read()
	while true do 
		local c = string.sub(input,pos,pos)
		if (c == "") then break end
		if (c == '"') then
			local text = ''
			local s,e,txt,c1
			repeat
				s,e,txt,c1 = string.find(input, '^(.-")(.)', pos+1)
				text = text..txt
				pos = e
				--print(txt, e, c1)
			until ((c1 == sep) or (c1 == "\r") or (c1 == "\n"))
			--print(string.sub(text,1,-2), c1)
			table.insert(row, string.sub(text,1,-2))
			c = c1
			pos = pos + 1
		else
			local s,e,text,c1 = string.find(input, "^([^%"..sep.."\r\n]-)([%"..sep.."\r\n])", pos)
			pos = e+1
			--print(text, c1)
			table.insert(row, text)
			c = c1
		end
		if c == "\n" then
			return row, pos
		end
		if c == "\r" then
			return row, pos+1
		end
	end
end

function slice_csv(csv, col)
  local sliced = {}

  for i = 1, #csv, 1 do
    sliced[#sliced+1] = csv[i][col]
  end

  return sliced
end

csv = {}
csv.load = function(filename, delimiter, header)
	local f,err = io.open(filename)
	if not f then
		print(err)
		return
	end
	local csv = f:read("*a")
	f:close()

	sep = string.sub(delimiter,1,1) or ','
	local pos = 1
	local t_csv = {}
	local f_header = nil
	local t_header = {}
	if header then
		t_header,pos = parse_row(csv, sep, pos)
		local head = {}
		for i,v in ipairs(t_header) do
			head[v] = i
		end
		f_header = function (t,k)
			local i = head[k]
			if i then
				return t[i]
			end
			return nil
		end
	end

	local row = {}
	row, pos = parse_row(csv, sep, pos)
	while row do
		if header then
			setmetatable(row, { __index = f_header })
		end
		table.insert(t_csv, row)
		row, pos = parse_row(csv, sep, pos)
	end
	return t_csv, t_header
end


return csv



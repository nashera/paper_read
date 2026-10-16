local bitexe = "bit.exe"
local bitpid = spawn(bitexe, "sync")
if not bitpid then
    print("bit.exe is not found")
    os.exit(1)
end
timeout = 10
capturelines = 1
local promptkeyword = "$"

while true do
    local rc = expect("Please provide a description of your changes [Enter 2 empty lines to finish]")
    local getTime = os.date("%c")
    if rc == 0 then
        sendln(getTime)
        sendln("\n")
        sleep(1)
        break
    else
        -- echo(string.format("Error keyword found \"$s\", Exit", _MATCH))
        break
    end
end

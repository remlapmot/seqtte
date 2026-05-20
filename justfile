# dev:
#     stata-mp "do developer"
# render:
#     quarto render site
#     rm -r docs/*
#     mv site/_site/* docs/
[working-directory: 'cscripts']
test:
    stata-mp -b "do master"
install:
    stata-mp -q <<< 'net install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace'
uninstall:
    stata-mp -q <<< 'ado uninstall seqtte'

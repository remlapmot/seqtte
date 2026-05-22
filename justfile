# dev:
#     stata-mp "do developer"
# render:
#     quarto render site
#     rm -r docs/*
#     mv site/_site/* docs/
[working-directory: 'cscripts']
test:
    stata-mp -b "do master"
in:
    printf 'ado uninstall seqtte\nnet install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace\nado dir seqtte\n' | stata-mp -q
un:
    printf 'ado uninstall seqtte\nado dir seqtte\n' | stata-mp -q

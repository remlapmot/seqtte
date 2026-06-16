# dev:
#     stata-mp "do developer"
# render:
#     quarto render site
#     rm -r docs/*
#     mv site/_site/* docs/

install: uninstall
    printf 'net install seqtte, from("https://raw.githubusercontent.com/remlapmot/seqtte/main/") replace\nado dir seqtte\n' | stata-mp -q

uninstall:
    #!/usr/bin/env bash
    # seqtte may be installed multiple times. `ado uninstall seqtte` fails with
    # r(111) ("matches more than one package") when there are duplicates, and each
    # uninstall renumbers the remaining packages, so we re-query `ado dir seqtte`
    # and remove one [#] at a time until none remain. The iteration cap is just a
    # safety stop against an unexpected non-terminating loop.
    for _ in $(seq 1 100); do
        n=$(printf 'ado dir seqtte\n' | stata-mp -q | grep -oE '\[[0-9]+\]' | head -1)
        [ -z "$n" ] && break
        printf 'ado uninstall %s\n' "$n" | stata-mp -q
    done
    echo "Remaining seqtte packages (should be none):"
    printf 'ado dir seqtte\n' | stata-mp -q

[working-directory('cscripts')]
test:
    printf 'cap noi ado uninstall seqtte\n' | stata-mp -q
    stata-mp -b "do master"

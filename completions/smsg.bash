# Bash completion for smsg (SAT Modulo Symmetries graph generator)

_smsg()
{
    local cur prev
    _init_completion || return

    case "$prev" in
        --dimacs|--simplify|--learned-clauses|--cube-file|--cube-file-test|--qcir|--sym-break-clauses|--lrat-output|--forbidden-subgraph-file|--forbidden-induced-subgraph-file)
            _filedir
            return
            ;;
        --vertices|-v|--frequency|--cutoff|--planarity-frequency|--frequency-forbidden-subgraphs|--min-chromatic-number|--chi|--coloring-algo|--triangle-vars|--max-learned-clause-size|--prerun|--assignment-cutoff|--simple-assignment-cutoff|--create-game-recursion-lvl|--cube-line|--cube-timeout|--timeout)
            return
            ;;
        --create-randomized-game)
            return
            ;;
        --initial-partition|--cubes-range|--cadical-config)
            return
            ;;
    esac

    if [[ "$cur" == -* ]]; then
        local opts="
            --help
            --all-graphs -a
            --vertices -v
            --directed -d
            --no-SMS
            --dimacs
            --hide-graphs
            --frequency
            --cutoff
            --initial-partition
            --colex-ordering
            --sym-break-clauses
            --connected -c
            --planar -p
            --planarity-frequency
            --forbidden-subgraph-file
            --forbidden-induced-subgraph-file
            --frequency-forbidden-subgraphs
            --min-chromatic-number --chi
            --coloring-algo
            --non-010-colorable
            --triangle-vars
            --simplify
            --learned-clauses
            --max-learned-clause-size
            --lrat-output
            --prerun
            --assignment-cutoff
            --simple-assignment-cutoff
            --create-game
            --create-randomized-game
            --create-game-recursion-lvl
            --cube-file
            --cube-file-test
            --cube-line
            --cubes-range
            --cube-timeout
            --timeout
            --cube-only-decisions
            --lookahead-only-edge-vars
            --cadical-config
            --qcir
            --polarity-hasing
        "
        COMPREPLY=($(compgen -W "$opts" -- "$cur"))
        return
    fi
} &&
complete -F _smsg smsg

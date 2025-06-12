package config      // CLI configuration file

// Tool + "flags" to allow specific callings
type Options struct {
    Tool   string
    Params map[string]string
}

func ParseArgs(args []string) Options {
    opts := Options{Params: make(map[string]string)}
    if len(args) > 0 {
        opts.Tool = args[0]
    }
    for _, arg := range args[1:] {
        kv := splitOption(arg)
        opts.Params[kv[0]] = kv[1]
    }
    return opts
}

// 
func splitOption(arg string) [2]string {
    var kv [2]string
    for i, ch := range arg {
        if ch == '=' {
            kv[0] = arg[:i]
            kv[1] = arg[i+1:]
            return kv
        }
    }
    kv[0] = arg
    kv[1] = ""
    return kv
}

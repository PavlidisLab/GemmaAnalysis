package ubic.gemma.script.example

import org.apache.commons.cli.*

class CliBuilder {
    String usage

    private Options options = new Options()
            .addOption('h', 'help', false, 'print this message')

    void usage() {
        def helpFormatter = new HelpFormatter();
        helpFormatter.printHelp(usage, options)
    }

    CommandLine parse(String[] strings) {
        try {
            def cl = new DefaultParser().parse(options, strings);
            if (cl.hasOption('h')) {
                usage()
                System.exit(0)
            }
            return cl
        } catch (ParseException e) {
            if (strings.contains("-h") || strings.contains("--help")) {
                usage()
                System.exit(0)
            }
            println e.message
            usage()
            System.exit(1)
        }
    }

    CliBuilder arg(Map args, String opt, String desc) {
        def builder = Option.builder(opt)
                .longOpt(args.longOpt as String)
                .argName(args.argName as String);
        if (args.required) {
            builder.required();
        }
        if (args.args) {
            builder.numberOfArgs(args.args as int)
        }
        if (args.required) {
            builder.required();
        }
        builder.desc(desc);
        options.addOption(builder.build())
        return this
    }
}

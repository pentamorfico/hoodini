import click


class MutuallyExclusiveOption(click.Option):
    """
    A custom click.Option that refuses to let certain other options be provided
    at the same time. Use it like:

        @click.option(
            "--foo",
            cls=MutuallyExclusiveOption,
            mutually_exclusive=["bar", "baz"],
            help="..."
        )
        @click.option(
            "--bar",
            cls=MutuallyExclusiveOption,
            mutually_exclusive=["foo"],
            help="..."
        )
    """

    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop("mutually_exclusive", []))
        super().__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.name in opts and opts.get(self.name) not in (None, False):
            for forbidden in self.mutually_exclusive:
                if forbidden in opts and opts.get(forbidden) not in (None, False):
                    msg = (
                        f"Illegal usage: '{self.name}' is mutually exclusive "
                        f"with '{forbidden}'. You cannot use both at once."
                    )
                    raise click.UsageError(msg)
        return super().handle_parse_result(ctx, opts, args)

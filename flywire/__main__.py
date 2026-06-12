"""Allow `python3 -m flywire ...` as an alternative entry point."""

from .cli import main

if __name__ == '__main__':
    main()

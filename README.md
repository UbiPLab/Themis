## File Integrations

This folder contains application integration files. They were tested on
Ubuntu 18.04 and 20.04.

For instructions how to build and run each application under Gramine, please
see the README.md in each subdirectory.

Please note that most of the files use oversimplified configurations which
are not secure. E.g., we frequently specify security-critical files as
`sgx.allowed_files`. If you take these files as templates for your own
production workloads, please inspect and harden the configurations.

We recommend to look at the (extensively commented) Redis file to get an idea
how to write the README, Makefile and manifest files. If you want to contribute
a new file to Gramine and you take the Redis file as a template, we
recommend to remove the comments from your copies as they only add noise (see
e.g. Memcached for a "stripped-down" file).

## Building Files

All our files use simple Makefiles to build the files and enable them
under Gramine. Use one of these commands:
- `make`: create non-SGX no-debug-log manifest
- `make DEBUG=1`: create non-SGX debug-log manifest
- `make SGX=1`: create SGX no-debug-log manifest
- `make SGX=1 DEBUG=1`: create SGX debug-log manifest

Use `make clean` to remove Gramine-generated artifacts and `make distclean` to
remove all build artifacts (if applicable).

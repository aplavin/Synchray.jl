# Contributing to Synchray.jl

Thank you for your interest in contributing to Synchray.jl!

## Reporting Bugs

If you encounter a bug, please open an issue on [GitHub Issues](https://github.com/aplavin/Synchray.jl/issues) with:
- A minimal reproducible example
- The Julia version and Synchray version you are using
- The expected vs. actual behavior

## Suggesting Features

Feature requests are welcome via GitHub Issues. Please describe the use case and, if possible, how you envision the API.

## Submitting Pull Requests

1. Fork the repository and create a feature branch.
2. Make your changes, keeping the existing code style (see below).
3. Add or update tests as appropriate.
4. Open a pull request with a clear description of the changes.

## Running Tests

As usual in Julia, type
```julia
] test
```
in the package environment.

To run a specific test item, run
```julia
@run_package_tests filter=ti -> ti.name == "test name"
```
in the test environment.

## Questions

For questions about usage or development, open a GitHub Issue or start a Discussion.

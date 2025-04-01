## Introduction

## Features


## Getting Started

### Prerequisites

ReadmeAI requires Python 3.9 or higher, and one of the following installation methods:

| Requirement                          | Details                          |
|--------------------------------------|----------------------------------|
| • [Python][python-link] ≥3.9         | Core runtime                     |
| **Installation Method** (choose one) |                                  |
| • [pip][pip-link]                    | Default Python package manager   |
| • [pipx][pipx-link]                  | Isolated environment installer   |
| • [uv][uv-link]                      | High-performance package manager |
| • [docker][docker-link]              | Containerized environment        |



### Usage


```sh
❯ readmeai --api openai -o readmeai-openai.md -r https://github.com/eli64s/readme-ai
```

> [!IMPORTANT]
> The default model set is `gpt-3.5-turbo`, offering the best balance between cost and performance.When using any model from the `gpt-4` series and up, please monitor your costs and usage to avoid unexpected charges.


```sh
❯ readmeai --repository /users/username/projects/myproject --api openai
```


### Testing

<!-- #### Using `pytest` [![pytest][pytest-shield]][pytest-link] -->

The [pytest][pytest-link] and [nox][nox-link] frameworks are used for development and testing.

Install the dependencies with uv:

```sh
❯ uv pip install --dev --group test --all-extras
```

Run the unit test suite using Pytest:

```sh
❯ make test
```

Using nox, test the app against Python versions `3.9`, `3.10`, `3.11`, and `3.12`:

```sh
❯ make test-nox
```



## Credits
Research Group:
- Adrian Batista-Planas
- Ernesto Quintas-Sanchez
- Richard Dawes (Advisor)

This work was partially supported by the Missouri University of Science and Technology’s Kummer Institute for Student Success and the United States Department of Energy (DOE), grant numbers DE-SC0019740 and DE-SC0025420.




# Distribution Instructions

## Pypi Upload instruction

- install `twine`
- `python setup.py sdist`
- Test environment
  - `python3 -m twine upload --repository testpypi dist/*`
  - `python3 -m pip install nemus-magpie --extra-index-url=https://test.pypi.org/simple/`
- Live Environment
  - `python3 -m twine upload dist/*` 
  - `python3 -m pip install nemus-magpie --extra-index-url=https://test.pypi.org/simple/`


## Documentation
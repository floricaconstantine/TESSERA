# Utility function to instantiate variables from elements in a named list. E.g., call this function on list(a=1, b=2) would result in variables a and b with values 1 and 2, respectively, in the calling frame/environment.

Utility function to instantiate variables from elements in a named list.
E.g., call this function on list(a=1, b=2) would result in variables a
and b with values 1 and 2, respectively, in the calling
frame/environment.

## Usage

``` r
variables_from_list(lst, target_environ = parent.frame())
```

## Arguments

- lst:

  List with named fields.

- target_environ:

  Environment to instantiate variables in.

## Value

Nothing; the global environment is modified.

## Note

This function is not currently used elsewhere in the codebase.

This function has identical functionality and essentially identical code
to extract.named in mvbutils.

## Author

Florica J Constantine, florica AT berkeley.edu

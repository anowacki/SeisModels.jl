# Contributing to SeisModels.jl

Thank you very much for your interest in supporting
SeisModels.jl with contributions.


## Contact

The primary maintainer of SeisModels.jl is Andy Nowacki
(@anowacki).  I welcome questions about the software,
which may be best raised as issues in the repo, but please
realise that I may not have much time to help you use the
package if your question is of the type: “How do I do …?”
I will do my best whenever I can.


## Code style guidelines and layout

If you open a pull request to the repo with fixes or new
functionality, I will do my best to help ensure it is
consistent in style with the rest of the package and then
merge it into the codebase.
However, it would be helpful if you can try and mirror
the style of the existing code as best you can to enable
your PR to be merged as quickly as possible.  This is best
done by looking at other parts of the package first.


## Dependencies

SeisModels.jl deliberately does not depend on many heavyweight
packages, to keep its loading and installation time small, and
therefore hopefully encourage experimentation with it.  Therefore,
please consider very carefully any new dependencies that you
may want to introduce to implement new functionality,
particularly any non-cross-platform ones.


## New functionality and models

Contributions of new models, model types, and functionality
to compute local and planetary properties are most welcome.
If you would like to add your own favourite model, please do
go ahead.

To suggest any changes to the code, please fork the repo,
clone it locally, make your changes, then open a pull request.


## Tests

Any new code (apart from documentation) **must** be accompanied
with tests to help make sure I don’t break your amazing work
when I fiddle with the package in future.


## Bugs

Reports of bugs—especially those leading to any inaccuracy
or incorrect values—are especially welcome.  To report one, please
[open an issue](https://github.com/anowacki/SeisModels.jl/issues/new)
at the repo and provide as much detail as possible to
reproduce the bug.  It is most helpful if you can provide
a [‘minimal (non-) working example’](https://stackoverflow.com/help/minimal-reproducible-example).


## Documentation improvements

PRs to improve the documentation are very welcome.  That includes
the online docs, docstrings and even this document.


## Licence

Please pay attention to the fact that the code is distributed
under the [MIT licence](https://en.wikipedia.org/wiki/MIT_License).
This means that this code can be freely re-used by anyone else
for any purpose so long as attribution is given to the authors.
Any code submitted to this repository will be licensed under the
MIT licence, hence you must be happy for this to be the case
for your contribution, and I will assume this is the case for
any code submitting to the repository.

It is also worth noting that the licence means one cannot simply
copy code that is licensed under something like the
[GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License)
unless the original authors agree to re-licence it as MIT for
this package.

## Coding style guide

It's important that Rivet code be fairly homogeneous in coding style, so that all the analyses and projections are good examples for non-core developers to copy, and so that everyone can read each class without having to double-think.

Some of the rules below might seem rather petty, but please stick by them (and complain at Andy if you think you've got a good reason for doing otherwise!) As with all rules, they're there to be broken, but only if you can justify doing so :-)

### General code style

 * Header files should have a `*.hh` suffix, implementation files with a `*.cc` suffix;
 * Private member variables should have a `_` prefix followed by a ''lower case'' letter, e.g. `private: int _myFoo;`. Don't use starting underscores followed by capital letters, or double underscores at all - both can clash with internal symbols;
 * Prefer the prefix incrementation to postfix incrementation: `++i` rather than `i++`. For e.g. iterators, this can be a not-insignificant performance boost, and since this is a well-known idiom it's lazy to do otherwise - when I see post-increments in code, it makes me wonder if the author really knows what they're doing!
 * Use the One True Brace Style (1TBS) for all code submitted to Rivet: see [wikipedia.](http://en.wikipedia.org/wiki/Indent_style#Variant:_1TBS) And also...
 * Your code will be more readable if you space things out a bit:
```
// prefer
if (foo < bar) {...
// to 
if(foo<bar){...
```
  , 
```
// prefer
for (int i = 0; i < N; ++i) {
// to 
for(int i=0;i<N;++i){
```
  and
```
// prefer
if (foo == bar) {
// to
if(foo==bar){
```

  It's not Fortran: spaces are not expensive! We will edit code not supplied in this format.
 * Also leave a couple of blank lines between function definitions, and occasional blank lines between blocks inside functions and class definitions to better emphasise the algorithm/data structures;
 * Use two spaces rather than tabs for indentation. You should be able to set up your code editor to do this automatically - for example, try adding this to your `.emacs` file if you're an emacs user:

```
(setq-default c-basic-indent 2)
(setq-default tab-width 4)
(setq-default indent-tabs-mode nil)
```

  While silly flame wars rage on this topic, most people are agreed that consistency is important: since most editors these days should be able to indent "to the correct amount" automatically, using spaces is completely unambiguous, with no extra typing required.

 * Although it may seem horrid to Europeans, use American spelling e.g. (color rather than colour) in function and publicly-visible parameter/variable names;
 * Use line wrapping sensibly rather than religiously - 80 characters is a good default length, but if wrapping your 110 character line (e.g. a for loop using lengthy iterator names) obscures the code structure then you'll need to make a choice. Maybe it can be reduced by using a sensible `typedef`, or a temporary variable - often sensible anyway;
 * Don't comment the `#endif`s in header guards - it just leads to inconsistencies when new header files are made by "copy and paste". Similarly, don't use comments which mention the filename: they lead to inconsistencies via renamings and copyings, and if you need to print the source then the IDE/editor will probably attach the filename to the printout anyway.
 * Try to minimise the {{{#include}}} entries in header files, especially if they refer to external library objects: binary library link dependencies are one thing, but header dependencies are usually one step too far. This may require storing member variables as pointers, so that a forward declaration of the variable's class can be used rather than a header `#include` that might induce a transitive dependency.
 * Only use the `auto` keyword where the true type is obscure, unknown, or too lengthy to type. When it is known that the type will be e.g. `Jet` or `Particle`, use that name explicitly: it makes the code far more readable.


### Analysis coding style

 * Prefer not to use pointers ''at all'' in projection and analysis code. The one exception to this is YODA histograms, where you have no option but to use pointers. See below for more information.
 * Histogram objects should be private and their names should start with "`_hist`" or "`_h`" (and similar for profiles and counters);
 * All Rivet analyses are plugins via the `Analysis` interface, and virtually no analyses inherit from any base class other than `Analysis`. Hence there is no need for header files, and analyses should be written completely inline with the implementations as part of the class definition, all in the `.cc` file.
 * Use the `DEFAULT_ANALYSIS_CTOR()` and similar macros for boilerplate code.
 * Analysis class names should be ALL-CAPS unless there is good reason not to.
 * Variables should start with lower-case letters, type names as MixedCase, and constants as ALL_CAPS.
 * Use the `const` keyword to protect your code against accidental modifications.
 * Pass complex objects by reference to any utility functions.
 * Loop over containers using the colon-separated "range for" loop, with a (const) reference as loop variable.
 * It is ok, perhaps even recommended, to make everything in an analysis code `public`, since there is no header, no client code, and therefore visibility levels have no effect.
 * Use the `MSG_DEBUG(...)`, `MSG_INFO(...)` etc. macros in place of `getLog() << ... << endl;`.
 * For the most part, do not use the `std::` prefix on STL types: they are already included without namespace into the `Analysis` scope.
 * Use the convenience types and functions to improve readability, e.g. `p.abseta()` rather than `fabs(p.eta())`, `Particles` and `Jets` rather than `vector<Particle>` and `vector<Jet>`, etc. (There are also `strings`, `doubles` and `ints` typedefs for `vector<those_types>`)
 * If your analysis includes the same sets of plots (binnings and cuts can differ), then don't register histograms for each energy, e.g. `_hist_blah_900GeV`, `_hist_blah_7000GeV`, etc.: just make one `_hist_blah` and use the `sqrtS()` function in the `init()` method of the analysis to book it from the appropriate histogram code. Then in the `analyze()` method, you can just call `fill()` without having to work out which variable you should be filling. This can save a ''lot'' of repetitive copy 'n' paste code, and we will reject supplied analyses which should do this and haven't, since otherwise they are a maintenance nightmare.


### Why to not use pointers

The "don't use pointers" rule may seem particularly perverse. Aren't they core to how C++ works? Well, for most purposes the answer is "no". While there are certainly areas of C++ code where pointers are useful, they tend to only be the places where references can't be used: polymorphic containers, storage of abstract base classes, and member variables with reference semantics. References are safe, while pointers are the single most common cause of segfaults and difficult-to-find bugs: prefer references whenever possible.

One particular reason to discourage pointers is that we want analysis and projection classes to be writeable by non-C++ experts: if pointers are involved, the level of required expertise is immediately raised (even if you don't realise that that's the case). Using pointers also forces projection authors to have to write custom constructors, destructors, copy constructors and copy assignment operators (cf. the "if you need one, you'll need all three" idiom): that's a lot of work and potential bugs that could have been avoided. This rule is sometimes referred to as the Law of the Big Three: see http://www.parashift.com/c++-faq-lite/coding-standards.html#faq-27.10 (and the rest of this excellent C++ resource!) for details.

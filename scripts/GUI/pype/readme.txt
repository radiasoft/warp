Readme/Help for PyPE (Python Programmer's Editor)
http://pype.sourceforge.net
http://come.to/josiah

PyPE is copyright (c) 2003 Josiah Carlson.

This software is licensed under the GPL (GNU General Public License) as it
appears here: http://www.gnu.org/copyleft/gpl.html
It is also included with this archive as gpl.txt.

The included STCStyleEditor.py, which is used to support styles, was released
under the wxWindows license and is copyright (c) 2001 - 2002 Riaan Booysen.
The copy included was also distributed with wxPython version 2.4.1.2 for
Python 2.2, and was not modified in any form.

The included stc-styles.rc.cfg was slightly modified from the original version
in order to not cause exceptions during style changes, and was also
distributed with wxPython version 2.4.1.2 for Python 2.2

The included wxProject.py was modified from the original version distributed
with wxPython in order to support being a drag and drop source for files.  I
do not know what license it was distributed under, but as the readme.txt says,
it was from a IBM developerWorks article at:
http://www-106.ibm.com/developerworks/library/l-wxpy/index.html
It does not seem to work when run with a version of wxPython with unicode
support.  I offer no support for wxProject, nor will I accept any bug reports
or patches to it.  I include it with PyPE because there was a message
requesting that a 'project manager' be included with PyPE.  So be it.  Now
there is one.  If someone wants to create a full-blown project manager out of
it, and distribute it separately, feel free, I won't mind, heck, I'll even
link you.

If you do not also receive a copy of gpl.txt with your version of this
software, please inform the me of the violation at the web page at the top of
this document.


Does anyone actually read this?  If so, please contact me, because I don't
think anyone actually does.  I like to be proven wrong.

#------------------------------- Requirements --------------------------------
PyPE has only been tested on Python 2.3 and wxPython 2.4.2.4.  It should work
on later versions of Python and wxPython...unless the namespace for wxPython
changes radically.

#-------------------------------- Visual bugs --------------------------------
The line and column indicator is accurate only after a keyboard key has been
released.  This is the result of EVT_STC_POSCHANGED not working properly.

#----------------------------------- Help ------------------------------------
How do I use the 'Todo' list?
1. Make sure that the 'Show Todo' has a check next to it in the 'Document'
   menu.
2. Restart PyPE to make sure that it has enabled the todo list.
3. Insert comments with any amount of whitespace preceeding them with the
   string 'todo:' in it (any capitalization of 'todo' is fine).  Everything to
   the right of 'todo:' will be the text in the todo list.  The number of
   exclamation points will be the '!' column.

The following lines are all valid todos;
# todo: fix the code below
        #todo:fix the code below!
    #        !TODo: fix the code below
# <I will not be noted in the todo message>tOdO:I will be

As you can see, it is pretty flexible.  As long as the line is a comment with
only whitespace before it and contains the string 'todo:' (with any
capitalization), you'll get an item in your todo list.



What is the difference between the fast and slow parser?
The fast one uses custom parser that basically splits the file into lines,
does a check to see if there is a function or class definition, then saves the
heirarchy information based on the level of indentation and what came before
it.  This can be innaccurate, as the fast parser will mistakenly believe that
the below function 'enumerate' is a method of MyException.  The slow parser
has no problems, it uses the standard Compile module.

class MyException(exceptions.Exception):
    pass
try:
    enumerate
except:
    def enumerate(inp):
        return zip(range(len(inp)), inp)

The fast parser also doesn't know anything about multi-line strings, so the
definition nada in the following line would be seen as a function, and not
part of a string.

old = 'this used to be a function\
def nada(inp):\
    return None'

Ah well, one has to give up something for speed.  Another thing given up is
that the fast parser will not pull out doc strings or handle multi-line
function definitions properly.  This may be changed in the future (parsing is
an entertaining and educational topic), but it will likely ALWAYS be
inaccurate.

Since 1.7, the slow parser is unused.  This is because the fast parser is
approximately 5-10 times faster and returns enough information to be useful.



Calltips:
How do you get usable Calltips?  Easy.  Hit F5.  This will also rebuild the
browsable source tree, autocomplete listing, and todo list (if it is enabled).



Autocompletion:
How do you get autocompletion?  Easy.  In the 'Document' menu, there is an
entry for 'Show autocomplete'.  Make sure there is a checkbox by it, and you
are set.  If you want to get a new listing of functions, hit the F5 key on
your keyboard.



Shell Commands:
I don't know how much other people use this feature, but I use it enough to
warrant the time I spent implementing it.  Basically this allows you to run
shell commands or scripts or even the script you are currently editing.  It
SHOULD just work for most things.

Play around with it, and remember to read the titles of the windows.



Code Snippets:
Ahh, what are code snippets?  Basically it is a saved-state multiple-entry
clipboard.  What is the use of that?  Well, let us say that you have a
template for interfaces to, let us say, commands for an interactive online
multiplayer game.  Each command needs to have a specific format, and with this
code snippet support, you don't need to switch to your template file, copy,
switch back and paste.  You can select your insertion point, and double click.
There are, of course, hot-keys for using code snippets while editing your
document.  Why?  Because I like having that option.  You can navigate and
insert code snippets without ever having your hands leave the keyboard.
Deleting a code snippet is a easy as making sure the listbox has keyboard
focus, and hitting 'delete' when the snippet you want to remove is selected.
Play around with it, you will (I believe) come to enjoy it.

Code Snippets are saved at program exit automatically.

If there is feedback/desire for being able to reorganize code snippets, I'll
add that support.



Bookmarked Paths:
Everyone will surely notice the menu for Pathmarks.  This menu allows you to
edit and access bookmarked paths with relative ease.  All it really does is
remember the paths you tell it to, and when you use one of the hotkeys or menu
items, it will change the current working directory to that new path. If you
attempt to open a file immediately afterwards, the open dialog will seek to
the path of the just used bookmark.  Nifty eh?  I like to think of it as being
able to have 'projects' without having to specify a project file.  I hate
project files.



CRLF/LF/CR line endings:
PyPE will attempt to figure out what kind of file was opened, it does this by
counting the number of different kind of line endings.  Which ever line ending
appears the most in an open file will set the line ending support for viewing
and editing in the window.  Also, any new lines will have that line ending.
New files will have the same line endings as the host operating system.

Additionally, copying from an open document will not change the line-endings.
Future versions of PyPE may support the automatic translation of text during
copying and pasting to the host operating system's native line endings.

Converting between line endings is now a menu item that is available in the
'Document' menu.



STCStyleEditor.py:
As I didn't write this, I can offer basically no support for it.  It seems to
work to edit python colorings, and if you edit some of the last 30 or so lines
of it, you can actually use the editor to edit some of the other styles that
are included.

If it just doesn't work for you, I suggest you revert to the copy of the
editor and stc-styles.rc.cfg that is included with the distribution of PyPE
you received.  As it is a known-good version, use it.



Expandable/collapseable/foldable code:
Since the beginning, there have been expandable and collapseale scopes thanks
to wxStyledTxtCtrl.  How to use them...
Given the below...
- class nada:
-     def funct(self):
-         if 1:
|             #do something
|             pass
Shift-clicking the '-' next to the class does this...
- class nada:
+     def funct(self):

Or really, it's like ctrl-clicking on each of the functions declared in the
scope of the definition.
Shift-clicking on the '-' a second time does nothing.
Shift-clicking on a '+' expands that item completely.

Control-clicking on a '+' or '-' collapses or expands the entirety of the
scopes contained within.

I don't know about you, but I'm a BIG fan of shift-clicking classes.  Yeah.
Play around with them, you'll get to loving how they work.



Find/Replace dialogs:
One big thing to note is how the find and find/replace dialogs work.  For
those of you who are annoyed with one's normal inability to enter in things
like newlines, this will be a great thing for you.

If you have ' or " as the first character in a find or find/replace dialog,
and what you entered is a proper string declaration in Python, that is the
string that will be found/replaced or replaced with.  Sorry, no raw string
support via r'c:\\new downloads\\readme.txt' so yeah.

As well, I've not implemented the ability to search up in the find dialog.
I may in future versions.  No matter what you have selected, it will always
search forward (and even wrap around if you keep hitting 'find next').

Converting between tabs and spaces:
So, you got tabs and you want spaces, or you have spaces and want to make them
tabs.  As it is not a menu option, you're probably wondering "how in the hell
am I going to do this".  Well, if you read the above stuff about the find and
replace dialogs, it would be trivial.
Both should INCLUDE the quotation marks.
To convert from tabs to 8 spaces per tab; replace "\t" with "        "
To convert from 8 spaces to one tab; replace "        " with "\t"

#-------------------------- How did PyPE come about --------------------------
The beginnings of PyPE was written from 10:30PM on the 2nd of July through
10:30PM on the 3rd of July.  Additional features were put together on the 4th
of July along with some bug fixing and more testing for version 1.0.
Truthfully, I've been using it to edit itself since the morning of the 3rd of
July, and believe it is pretty much feature-complete (in terms of standard
Python source editing).  There are a few more things I think it would be nice
to have, and they will be added in good time.

On the most part, this piece of software should work exactly the way you
expect it to.  That is the way I wrote it.  As a result, you don't get much
help in using it (mostly because I am lazy).  When questions are asked, I'll
add the question and answer into the FAQ, which is at the end of this
document.


The majority of the things that this editor can do are in the menus.  Hot-keys
for things that have them are listed next to their menu items.  As I am still
learning all the neat things one can do with wxStyledTxtCtrl, I don't know all
the built-in features, and this is likely as much of a learning experience for
me as you.


#------------------------------------ FAQ ------------------------------------
When I use 'Shell'->'Run current file', my script doesn't run.  WTF?
If you were to open a dos/shell prompt in the same path as your script, then
type in the name of your script, would it run?  PyPE will only run scripts
that are able to be run at such a prompt.

The most common of these that seem like they should run, but don't, is when a
file name or path has unicode characters in them.  os.system and os.spawnv
don't handle unicode very well.  os.system crashes.  os.spawnv automatically
translates unicode characters into question mark '?' characters.  Hell, even
os.listdir converts unicode characters into '?'.



What's the deal with the version numbering scheme?
Early in development, PyPE raised version numbers very quickly.  From 1.0 to
1.5, not much more than 2 months passed.  In that time, most of the major
architectural changes that were to happen, happened.  This is not the reason
for the version number change.  Really it was so that the MAJOR versions could
have their own point release (1.0 being the first), and minor bugfixes on the
point releases would get a minor release number (like 1.0.1).

Then, at around PyPE 1.4.2, I had this spiffy idea.  What if I were to release
a series of PyPE versions with the same version numbers as classic Doom?  I
remembered updating to 1.1, then to 1.2a, etc.  My favorite was 1.666.  Ah hah!
PyPE 1.6.6.6, the best version of PyPE ever.

I decided that I would slow version number advancement, if only so that people
didn't get sick of new releases of PyPE being numbered so much higher, being
that there are other projects where version 1.0 (if it is ever released) would
be the perfect release, that there is no more changes to be made.

Then the more I thought about it, the more I realized that it doesn't matter
at all, I mean, Emacs is on version 20+.  *shrug*

Yeah.  To be honest, version numbers advance as my whims decree.  As of this
writing, PyPE 1.6.5.2 was released yesterday.  1.6.5.2 has a nice ring to it.
Maybe the next version will be 1.6.5.3, or 1.6.6, or 1.7.  Who knows?  It
doesn't really matter to tell the truth, as long as bugs get fixed, features
are added, and people keep using it, what does it matter what version you are
using?

#-------------------------------- Thank Yous ---------------------------------
Certainly there are some people I should thank, because without them, the
piece of software you are using right now, just wouldn't be possible.

Guido van Rossum - without Guido, not only would I not have Python, I also
wouldn't have had some of the great inspiration that IDLE has offered.  IDLE
is a great editor, has some excellent ideas in terms of functionality, but it
unfortunately does not offer the extended functionality I want, and it hurts
my brain to use tk, so I cannot add it myself.  Guido, my hat goes off to you.

The people writing wxWindows and wxPython - without you, this also would not
have been possible.  You have made the most self-consistent GUI libraries that
I have ever used, made them easy to use, and offer them on every platform that
I would ever want or need.  You rock.

The people writing Scintilla - as wxStyledTextCtrl is a binding for scitilla
for wxWindows, which then has bindings for wxPython, basically ALL the REAL
functionality of the editor you are now using is the result of Scintilla.  The
additional things like tabbed editing, hotkeys, etc., they are mere surface
decorations in comparison to what it would take to write everything required
for a text editor from scratch.  Gah, an editor widget that just works?  Who
would have figured?

To everyone who I have already thanked: thank you for making PyPE an almost
trivial task.  It would have been impossible to go so far so fast by hand in
any other language using any other GUI toolkit or bindings.

My wife - because without her, I would likely be a pathetic shell of a man.
#------------------------------- End of file. --------------------------------

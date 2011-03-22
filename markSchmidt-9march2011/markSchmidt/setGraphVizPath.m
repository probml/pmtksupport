PATH = getenv('PATH');
PATH = strcat(PATH,':/usr/local/bin');
setenv('PATH',PATH);
getenv('PATH')
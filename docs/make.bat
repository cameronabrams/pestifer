@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=source
set BUILDDIR=build
set APIDOC_OUT=%SOURCEDIR%\api
set APIDOC_SRC=..\pestifer

if "%1" == "" goto help
if "%1" == "apidoc" goto apidoc
if "%1" == "config-ref" goto config_ref
if "%1" == "html" goto html
if "%1" == "clean" goto clean

%SPHINXBUILD% >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
	echo.installed, then set the SPHINXBUILD environment variable to point
	echo.to the full path of the 'sphinx-build' executable. Alternatively you
	echo.may add the Sphinx directory to PATH.
	echo.
	echo.If you don't have Sphinx installed, grab it from
	echo.http://sphinx-doc.org/
	exit /b 1
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:apidoc
if exist %APIDOC_OUT% rmdir /s /q %APIDOC_OUT%
sphinx-apidoc -f -M -e --implicit-namespaces --tocfile API -o %APIDOC_OUT% %APIDOC_SRC%
goto end

:config_ref
cd %SOURCEDIR%
yclept make-doc ..\..\pestifer\schema\base.yaml --root config_ref --footer-style raw-html
cd ..
goto end

:html
if exist %APIDOC_OUT% rmdir /s /q %APIDOC_OUT%
sphinx-apidoc -f -M -e --implicit-namespaces --tocfile API -o %APIDOC_OUT% %APIDOC_SRC%
%SPHINXBUILD% -M html %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
goto end

:clean
if exist %BUILDDIR% rmdir /s /q %BUILDDIR%
if exist %APIDOC_OUT% rmdir /s /q %APIDOC_OUT%
goto end

:help
%SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%

:end
popd

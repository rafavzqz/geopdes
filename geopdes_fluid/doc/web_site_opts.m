%% This is the script used to produce the online function reference web pages
%% that you can browse at http://geopdes.sourceforge.net/functions/geopdes_fluid
%% the script is included in the distribution in case you may need a local copy
%% to use when the internet is unreachable.
%% To generate html docs:
%% 1 - install and load the geopdes_fluid and generate_html packages
%% 2 - cd to this directory and run this script from within octave

ws_version = "forge";
site_root  = "http://geopdes.sf.net";
pkg_name   = "geopdes_fluid"

if strcmpi (ws_version, "docbrowser")

opt.header ="\n\
<!DOCTYPE html PUBLIC ""-//W3C//DTD XHTML 1.0 Strict//EN""\n\
 ""http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"">\n\
<html xmlns=""http://www.w3.org/1999/xhtml"" lang=""en"" xml:lang=""en"">\n\
  <head>\n\
  <meta http-equiv=""content-type"" content=""text/html; charset=iso-8859-1"" />\n\
  <meta name=""date"" content=""26-May-2010""/>\n\
  <meta name=""author"" content=""The Octave Community"" />\n\
  <title>%title</title>\n\
  <link rel=""stylesheet"" type=""text/css"" href=""doc.css"" />\n\
  </head>\n\
<body>\n\
<div id=""top"">Function Reference</div>\n\
<div id=""doccontent"">";

opt.footer="\n\
</div>\n\
</body>\n\
</html>";

opt.include_demos = true;
opt.include_overview = true;

opt.overview_title="Overview: %name";
opt.title="Function: %name";

elseif strcmpi (ws_version, "forge")

opt = get_html_options ("octave-forge");

opt.include_package_news = 0;

%% opt.download_link = "http://downloads.sourceforge.net/geopdes/%name-%version.tar.gz?download";
opt.download_link = "http://sourceforge.net/project/geopdes/files/GeoPDEs_full.tar.gz/download";

opt.header = sprintf("<!DOCTYPE html PUBLIC ""-//W3C//DTD XHTML 1.0 Strict//EN""\n\
 ""http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"">\n\
<html xmlns=""http://www.w3.org/1999/xhtml"" lang=""en"" xml:lang=""en"">\n\
  <head>\n\
  <meta http-equiv=""content-type"" content=""text/html; charset=iso-8859-1"" />\n\
  <meta name=""date"" content=""%s""/>\n\
  <meta name=""author"" content=""%s"" />\n\
  <title>%%title</title>\n\
  <link rel=""stylesheet"" type=""text/css"" href=""%s/functions/docs.css"" />\n\
  </head>\n\
<body>\n\
<div id=""top""><a href=""http://geopdes.sourceforge.net"">Home</a> &middot; <a href=""%s/functions/%s/overview.html"">Function Reference</a> </div>\n\
<div id=""doccontent"">", date, "rafavzqz@users.sourceforge.net", site_root, site_root, pkg_name);

opt.footer= sprintf("\n\
<div id=""left_footer"">\n\
  page maintained by <address><a href=""mailto:rafavzqz@users.sourceforge.net""> rafavzqz@users.sourceforge.net </a></address>\n\
  <!-- hhmts start --> Last modified: %s <!-- hhmts end -->\n\
</div>\n\
<div id=""right_footer"">\n\
<a href=""http://sourceforge.net/projects/geopdes""><img src=""http://sflogo.sourceforge.net/sflogo.php?group_id=352940&amp;type=9"" width=""80"" height=""15"" alt=""Get geopdes at SourceForge.net. Fast, secure and Free Open Source software downloads"" /></a> </div>\n\
\n\
</body>\n\
</html>", date);
endif

generate_package_html (pkg_name, [pkg_name "_html"], opt);

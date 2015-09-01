#!/usr/bin/perl
use strict;
use warnings;

use LWP;

my $browser = LWP::UserAgent->new;
my $url = 'http://www.acme.com/';

my $response = $browser->get($url, 'User-Agent' => 'Mozilla/4.0 (compatible; MSIE 7.0)',);
die "Error getting $url" unless $response->is_success;
print 'Content type is ', $response->content_type;
print 'Content is:';
print $response->content;


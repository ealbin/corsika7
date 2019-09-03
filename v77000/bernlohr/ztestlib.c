/*    this test function is used by ./configure to test whether this flag 
 *    was used or not.
*/
/*
 *    If iact.c is used but atmo.c is not, then compile iact.c with
 *    the -DNO_EXTERNAL_ATMOSPHERES compiler switch.
*/
#ifndef NO_EXTERNAL_ATMOSPHERES
char testatmo()
{ return 0;
}
#else
char testnoatmo()
{ return 0;
}
#endif
/*
 *    test the -DIACTEXT compiler switch.
*/
#ifdef IACTEXT
char testext()
{ return 0;
}
#else
char testnoext()
{ return 0;
}
#endif

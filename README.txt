*******************************************************************************
The development model
*******************************************************************************

The cmiv-dira development model uses the "integration-manager" workflow
(Chacon, 2009, http://git-scm.com/book/en/Distributed-Git-Distributed-Workflows#Integration-Manager-Workflow):

1. The project maintainer pushes to their public repository.
2. A contributor clones that repository and makes changes.
3. The contributor pushes to their own public copy.
4. The contributor sends the maintainer an e-mail asking them to pull changes.
5. The maintainer adds the contributor’s repo as a remote and merges locally.
6. The maintainer pushes merged changes to the main repository.

The Git repository on code.google.com uses two long-running branches:
 - master
 - develop
(Chacon, 2009, http://git-scm.com/book/en/Git-Branching-Branching-Workflows#Long-Running-Branches)

*******************************************************************************
References:
*******************************************************************************
Chacon Scott, Pro Git (Expert's Voice in Software Development), accessed:
2014-07-10

##################
The Main Page
##################

This is the main page on which the intro to the documentation lives.

******************
A heading
******************
Some text to fill in a gap

.. tip::
    A tip for you

.. code-block:: Fortran

    subroutine test(field)
        type(field), intent(inout) :: field
        
        field = 2
    end subroutine
    
.. image:: images/fieldordering.svg
    :width: 800
    

.. sidebar:: Topic Title
    
    Some information about a topic.
    
.. math::

    \partial/\partial t = \nabla x

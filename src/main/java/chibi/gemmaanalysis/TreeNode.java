/*
 * The Gemma project
 *
 * Copyright (c) 2007 Columbia University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
package chibi.gemmaanalysis;

import java.io.Serializable;

import cern.colt.list.ObjectArrayList;

/**
 * TODO Document Me
 *
 * @author Paul
 * @version $Id: TreeNode.java,v 1.3 2013/09/26 22:13:27 paul Exp $
 */
public class TreeNode implements Comparable<TreeNode>, Serializable {

    private static final long serialVersionUID = 7454412470040843963L;
    public static final int MASKBITS = 1;
    public static final int COMMONBITS = 2;
    public static final int LEVEL = 3;
    public static final int ORDER = 4;
    private static int SORTING = MASKBITS;

    /**
     *
     */
    public static void reSetSorting() {
        TreeNode.SORTING = TreeNode.MASKBITS;
    }

    /**
     * @param sorting
     */
    public static void setSorting( int sorting ) {
        TreeNode.SORTING = sorting;
    }

    private long id;
    private Integer maskBits = 0;
    private long[] mask;
    private TreeNode closestNode = null;
    private Integer commonBits = 0;
    private ObjectArrayList children;
    private TreeNode parent = null;

    private Integer level = 0;

    private Integer order = 0; // for tree generation

    /**
     * @param id
     * @param mask
     * @param child
     */
    public TreeNode( long id, long[] mask, ObjectArrayList child ) {
        this.id = id;
        this.children = child;
        this.mask = mask;
        this.maskBits = LinkMatrix.countBits( this.mask );
    }

    /*
     * (non-Javadoc)
     *
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    @Override
    public int compareTo( TreeNode o ) {
        int res = 0;
        switch ( TreeNode.SORTING ) {
            case TreeNode.MASKBITS:
                res = maskBits.compareTo( o.maskBits ) * ( -1 );
                break;
            case TreeNode.COMMONBITS:
                res = commonBits.compareTo( o.commonBits ) * ( -1 );
                break;
            case TreeNode.LEVEL:
                res = level.compareTo( o.level ) * ( -1 );
                break;
            case TreeNode.ORDER:
                res = order.compareTo( o.order ) * ( -1 );
                break;
            default:
                break;
        }
        return res;
        // return maskBits.compareTo(o.maskBits)*(-1);
    }

    public ObjectArrayList getChildren() {
        return children;
    }

    public TreeNode getClosestNode() {
        return closestNode;
    }

    public Integer getCommonBits() {
        return commonBits;
    }

    public long getId() {
        return id;
    }

    public Integer getLevel() {
        return level;
    }

    public long[] getMask() {
        return mask;
    }

    public Integer getMaskBits() {
        return maskBits;
    }

    public Integer getOrder() {
        return order;
    }

    public TreeNode getParent() {
        return parent;
    }

    /**
     * @param child
     */
    public void setChildren( ObjectArrayList child ) {
        this.children = child;
    }

    /**
     * @param closestNode
     */
    public void setClosestNode( TreeNode closestNode ) {
        this.closestNode = closestNode;
        commonBits = LinkMatrix.overlapBits( mask, closestNode.mask );
    }

    /**
     * @param level
     */
    public void setLevel( int level ) {
        this.level = level;
    }

    /**
     * @param allBits
     */
    public void setMask( long[] allBits ) {
        this.mask = allBits;
    }

    /**
     * @param order
     */
    public void setOrder( int order ) {
        this.order = order;
    }

    /**
     * @param parent
     */
    public void setParent( TreeNode parent ) {
        this.parent = parent;
    }

}

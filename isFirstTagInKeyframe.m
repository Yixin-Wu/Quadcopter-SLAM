function a = isFirstTagInKeyframe(tag_id, tag_and_keyframe)
keyframe = tag_and_keyframe(tag_and_keyframe(:,1) == tag_id,2);
tag_id_in_this_keyframe = tag_and_keyframe(tag_and_keyframe(:,2) == keyframe, 1);
if tag_id == tag_id_in_this_keyframe(1)
    a = 1;
else
    a = 0;
end
end